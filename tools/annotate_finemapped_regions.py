# to annotate finemapped regions back to window's annotations
import pandas as pd
import deepdish as dd
from pybedtools import BedTool
import sys


if __name__ == '__main__':
    finemapped_sites = BedTool(sys.argv[1])
    finemapped_sites_df = finemapped_sites.to_dataframe()
    
    ranking = pd.read_csv(sys.argv[2], sep = '\t')

    # read annotation
    window = pd.read_csv(sys.argv[3],
                        sep = '\t')
    window_bed = BedTool.from_dataframe(window)

    outf = sys.argv[4]

    # find windows overlapping finemapped regions
    site_annotation_finemap = finemapped_sites.intersect(
        window_bed,s = True, wb = True).to_dataframe(
        names = finemapped_sites_df.columns.tolist()+[
            i+'.1' for i in window.columns.tolist()[:6]]+window.columns.tolist()[6:])
    
    site_annotation_finemap['overlap_size'] = site_annotation_finemap['end']-site_annotation_finemap['start']
    site_annotation_finemap['gene_types_rank'] = site_annotation_finemap['gene_types'].map(ranking.set_index('accession')['rank'])

    # get unique 1-to-1 mapping, priorize by ranking and overlap size (middle of peak)
    site_annotation_uniq = site_annotation_finemap.sort_values(
        by = 'overlap_size', ascending = False
                                   ).sort_values(by = 'gene_types_rank', ascending = True
                                   ).drop_duplicates(subset = 'name', keep = 'first')
    
    # save to file
    site_annotation_uniq.to_csv(outf, sep = '\t', index = False)