# ensure gcc is installed 
# create rskipper conda environment from rskipper.yml
# conda env create -f rskipper.yml
# set paths
CONDA_DIR=$1 # CONDA_DIR=/home/eboyle/miniconda3
R_DIR=$2 # R_DIR=/projects/ps-yeolab3/eboyle/encode/pipeline/gran

# load conda environment for R dependencies
conda activate rskipper

# make directories for downloading source
mkdir -p ${R_DIR}/src
mkdir ${R_DIR}/texlive

# set library files based on conda directory
export LD_LIBRARY_PATH=${CONDA_DIR}/envs/rskipper/lib/:$LD_LIBRARY_PATH

cd ${R_DIR}/src # working directory of your choice
wget --no-check-certificate https://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz # or curl instead of wget
zcat < install-tl-unx.tar.gz | tar xf -
cd install-tl-*
export TEXLIVE_INSTALL_PREFIX=${R_DIR}/texlive
perl ./install-tl --no-interaction # as root or with writable destination
# Finally, prepend /usr/local/texlive/YYYY/bin/PLATFORM to your PATH,
export PATH=${R_DIR}/texlive:$PATH # e.g., /usr/local/texlive/2023/bin/x86_64-linux

# install curl into conda directory manually (version on conda forge does not work)
cd ${R_DIR}/src
git clone https://github.com/curl/curl.git
cd curl
git checkout tags/curl-7_86_0
./buildconf
./configure --with-openssl --prefix=${CONDA_DIR}/envs/rskipper/
make
make install
cd ..

# install base R
curl https://cran.r-project.org/src/base/R-4/R-4.1.3.tar.gz > R-4.1.3.tar.gz
tar -xf R-4.1.3.tar.gz
cd R-4.1.3
./configure --enable-R-shlib --enable-shared --with-x=no --prefix=${R_DIR} CPPFLAGS="-I${CONDA_DIR}/envs/rskipper/include/" LDFLAGS="-L${CONDA_DIR}/envs/rskipper/lib/" && \
make && \
make install

# install versioned R packages using groundhog and BiocManager
cd ${R_DIR}/bin
./R -e "dir.create(Sys.getenv('R_LIBS_USER'),recursive=TRUE)"
./R -e "install.packages('groundhog',repos = 'http://cran.us.r-project.org')"
./R -e "groundhog::groundhog.library(c('tidyverse', 'VGAM', 'viridis', 'ggrepel', 'RColorBrewer', 'Rtsne', 'ggupset', 'ggdendro', 'cowplot', 'BiocManager'), '2022-03-11', force.install=TRUE)" && ./R -q -e "BiocManager::install(c('GenomicRanges','fgsea','rtracklayer'))"
