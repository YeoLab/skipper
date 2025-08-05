from dotenv import load_dotenv
import os
from huggingface_hub import login
# Imports
from transformers import AutoTokenizer, TrainingArguments, Trainer, AutoModelForSequenceClassification
import torch
import matplotlib.pyplot as plt
import numpy as np
from accelerate.test_utils.testing import get_backend
import warnings
from peft import LoraConfig, TaskType
from peft import get_peft_model
from pathlib import Path
import seqdata as sd
from datasets import Dataset
from scipy.stats import pearsonr
import sys
import torch

def to_dict(sdata):
    """
    Convert sdata to compatible format

    Parameters
    ----------
    test_sdata
        seqdata object with attribute seq and name
    output_filepath : str
        Path to write the FASTA file to.
    """
    seqs = []
    for i in range(len(sdata["name"])):
        s = "".join([c.decode() for c in sdata["seq"][i].values.tolist()])

        seqs.append(s)
    # Write all records in one go
    pseudocount = 10
    window_odd = (sdata["n_IP"] + pseudocount * sdata["gc_fraction"]) / (
        sdata["n_IN"] + pseudocount * (1 - sdata["gc_fraction"])
    )
    gc_bin_odd = sdata["gc_fraction"] / (1 - sdata["gc_fraction"])
    l2or = np.log2(window_odd) - np.log2(gc_bin_odd)

    records = {}
    records["data"] = seqs
    records["labels"] = l2or.values
    return records

def tokenize_function(examples):
    """Tokenize the input sequences using the tokenizer."""
    outputs = tokenizer(examples["data"])
    return outputs

def compute_metrics_pearson(eval_pred):
    """Compute the Pearson correlation coefficient between predictions and labels."""
    references = eval_pred.label_ids
    r, _ = pearsonr(eval_pred.predictions.squeeze(-1), eval_pred.label_ids)
    return {"d_log_odds_pearson": r}

if __name__ == "__main__":

    skipper_path = Path(sys.argv[1])
    exp = sys.argv[2]
    model_name = sys.argv[3] #"nucleotide-transformer-500m-human-ref"
    # Load environment variables from .env from working directory
    # go to huggingface and create a token and save it in .env
    load_dotenv()

    # Get the token
    hf_token = os.getenv("HF_TOKEN")

    # Login using the token
    login(token=hf_token)
    device, _, _ = get_backend()

    # set num labels to 1 to trigger regression

    # Load the model
    model = AutoModelForSequenceClassification.from_pretrained(
        f"InstaDeepAI/{model_name}", num_labels=1
    )
    model = model.to(device)
    if device != 'cuda':
        warnings.warn('You are not using a GPU, this will be very slow and never finish haha!')

    # set up LoRA finetuning
    peft_config = LoraConfig(
        task_type=TaskType.FEATURE_EXTRACTION,
        inference_mode=False,
        r=1,
        lora_alpha=32,
        lora_dropout=0.1,
        target_modules=["query", "value"],
        # modules_to_save=["intermediate"] # modules that are not frozen and updated during the training
    )

    lora_classifier = get_peft_model(
        model, peft_config
    )  # transform our classifier into a peft model
    lora_classifier.print_trainable_parameters()
    lora_classifier.to(device)  # Put the model on the GPU

    # Load data prepared for RBPNet
    data_dir = skipper_path / f"output/ml/rbpnet_data/{exp}"
    train_sdata = sd.open_zarr(data_dir / "train.zarr").load()
    valid_sdata = sd.open_zarr(data_dir / "valid.zarr").load()
    test_sdata = sd.open_zarr(data_dir / "test.zarr").load()


    # sys.path.append("/tscc/nfs/home/hsher/projects/RBPNet/")
    ds_train_promoter = Dataset.from_dict(to_dict(train_sdata))
    ds_validation_promoter = Dataset.from_dict(to_dict(valid_sdata))
    ds_test_promoter = Dataset.from_dict(to_dict(test_sdata))

    # Load the tokenizer
    tokenizer = AutoTokenizer.from_pretrained(
        f"InstaDeepAI/{model_name}"
    )
    # # Creating tokenized promoter dataset
    tokenized_datasets_train_promoter = ds_train_promoter.map(
        tokenize_function,
        batched=True,
        remove_columns=["data"],
    )
    tokenized_datasets_validation_promoter = ds_validation_promoter.map(
        tokenize_function,
        batched=True,
        remove_columns=["data"],
    )
    tokenized_datasets_test_promoter = ds_test_promoter.map(
        tokenize_function,
        batched=True,
        remove_columns=["data"],
    )

    # Path to save the model
    save_path = (
        skipper_path
        / f"output/ml/nt_lora/{exp}/{model_name}"
    )
    save_path.mkdir(exist_ok=True, parents=True)
    batch_size = 8
    args_promoter = TrainingArguments(
        save_path,
        remove_unused_columns=False,
        eval_strategy="steps",
        save_strategy="steps",
        learning_rate=5e-4,
        per_device_train_batch_size=batch_size,
        gradient_accumulation_steps=1,
        per_device_eval_batch_size=64,
        num_train_epochs=2,
        logging_steps=100,
        load_best_model_at_end=True,  # Keep the best model according to the evaluation
        metric_for_best_model="d_log_odds_pearson",
        label_names=["labels"],
        dataloader_drop_last=True,
        max_steps=2000,
    )

    trainer = Trainer(
        # model.to(device),
        lora_classifier,
        args_promoter,
        train_dataset=tokenized_datasets_train_promoter,  # to fix
        eval_dataset=tokenized_datasets_validation_promoter,
        tokenizer=tokenizer,
        compute_metrics=compute_metrics_pearson,
    )

    train_results = trainer.train()
    lora_classifier.save_pretrained(save_path / "final_state")

    # plot training curve
    curve_evaluation_f1_score = [
        [a["step"], a["eval_d_log_odds_pearson"]]
        for a in trainer.state.log_history
        if "eval_d_log_odds_pearson" in a.keys()
    ]
    eval_f1_score = [c[1] for c in curve_evaluation_f1_score]
    steps = [c[0] for c in curve_evaluation_f1_score]
    plt.plot(steps, eval_f1_score, "b", label="pearson correlation (d_log_odds)")
    plt.title(f"Validation pearson correlation for {exp} prediction")
    plt.xlabel("Number of training steps performed")
    plt.ylabel("Validation pearson correlation score")
    plt.legend()
    plt.show()
    plt.savefig(save_path / "training_curve.png")

    

    torch.cuda.empty_cache()
    all_eval = {}
    test_score = trainer.predict(tokenized_datasets_test_promoter).metrics[
        "test_d_log_odds_pearson"
    ]
    all_eval["all"] = test_score

    import pandas as pd

    coverage = pd.DataFrame(
        {
            "n_IP": torch.from_numpy(test_sdata["n_IP"].values),
            "n_IN": torch.from_numpy(test_sdata["n_IN"].values),
            "name": torch.from_numpy(test_sdata["name"].values),
        }
    )
    coverage["total"] = coverage["n_IP"] + coverage["n_IN"]
    # get genomic bin  for each coverage threshold
    nread2sdata = {}
    for nread in [10, 50, 100, 200]:
        nread2sdata[nread] = test_sdata.sel(
            _sequence=((test_sdata["n_IP"] + test_sdata["n_IN"]) > nread)
        )
    nread2tokenized = {}
    for nread in nread2sdata:
        sdata = nread2sdata[nread]
        dataset_sub = Dataset.from_dict(to_dict(sdata))

        tokenized_data_sub = dataset_sub.map(
            tokenize_function,
            batched=True,
            remove_columns=["data"],
        )
        dlogcorr = trainer.predict(tokenized_data_sub).metrics["test_d_log_odds_pearson"]
        all_eval[nread] = dlogcorr
        print(nread, dlogcorr)
    pd.Series(all_eval).to_csv(save_path / "d_log_odds_corr.csv")