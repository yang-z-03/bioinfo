zsh theme: skaro

conda init zsh

vim .condarc:
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  
conda install mamba -c conda-forge

mamba create -n omicverse.py310.r43 python=3.10 r-base=4.3 pkg-config
eval '$(mamba shell hook --shell zsh)'
mamba shell init --shell zsh --root-prefix=~/miniconda
mamba activate omicverse.py310.r43
vim ~/miniconda/envs/omicverse.py310.r43/conda-meta/pinned

# =====================================================================================

# install a user version of gcc compiler
# mamba install gcc -c conda-forge
# torch distributions should be managed using mamba.
mamba install pytorch torchvision torchaudio cpuonly -c pytorch
mamba install pyg -c pyg
# omicverse is managed using pip.
python -m pip install -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple --upgrade pip
pip config set global.index-url https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
pip install -U omicverse
pip install -U numba
pip install scvi-tools
# integrate into jupyter lab.
pip install ipykernel
python -m ipykernel install --user --name 'omicverse.python' --display-name 'omicverse.py310'
# installed kernelspec omicverse.python in 
# /home/data/yangzhen/.local/share/jupyter/kernels/omicverse.python

mkdir /home/data/yangzhen/.config/r-specs
mkdir /home/data/yangzhen/.config/r-specs/omicverse.r
touch /home/data/yangzhen/.config/r-specs/omicverse.r/.profile

# required by r package `gert`
mamba install libgit2
# required by r package `tidyverse`
mamba install libiconv 

# configure R.
R:
installed.packages() |> rownames()
options('repos' = c(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN'))
install.packages('IRkernel')
IRkernel::installspec(
  name = 'omicverse.r.torch-260',
  displayname = 'omicverse.r43 (torch-2-6-0)',
  verbose = TRUE
)
# installed kernelspec omicverse.r in 
# /home/data/yangzhen/.local/share/jupyter/kernels/omicverse.r

# for r 4.3, you should install matrix 1.6-5.
install.packages('~/shared/r-pkg/matrix-1-6-5.tar.gz', repos = NULL, Ncpus = 100, type = 'source')
install.packages('~/shared/r-pkg/mass-7-3-60.tar.gz', repos = NULL, Ncpus = 100, type = 'source')
# bioconductor version of fgsea contain a dependency problem. use the development branch
install.packages('tidyverse', Ncpus = 100)
BiocManager::install(c('devtools', 'remotes'), ask = FALSE, update = FALSE, Ncpus = 100)
install_github("ctlab/fgsea") 

install.packages('BiocManager', ask = FALSE, update = FALSE)
options(BioC_mirror = 'https://mirrors.westlake.edu.cn/bioconductor')
BiocManager::install('3.20')
required = c('Seurat', 'SingleR', 'clusterProfiler', 'crayon', 'tibble')
BiocManager::install(required, ask = FALSE, update = FALSE, Ncpus = 100)

options('repos' = c(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN'))
options(BioC_mirror = 'https://mirrors.westlake.edu.cn/bioconductor')


# =====================================================================================
# ============================ pip pytorch environment ================================
# =====================================================================================

mamba create -n omicverse.py310.r43.torch113 python=3.10 r-base=4.3 pkg-config
mamba activate omicverse.py310.r43.torch113
vim ~/miniconda/envs/omicverse.py310.r43.torch113/conda-meta/pinned
# pip install torch==1.13.1+cpu torchvision==0.14.1+cpu torchaudio==0.13.1 --extra-index-url https://download.pytorch.org/whl/cpu
pip config set global.index-url https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple
pip install -U omicverse
pip install -U numba
pip install scvi-tools
# integrate into jupyter lab.
pip install ipykernel
python -m ipykernel install --user --name 'omicverse.python.torch-260' --display-name 'omicverse.py310 (torch-2-6-0)'

# required by r package `gert`
mamba install libgit2
# required by r package `tidyverse`
mamba install libiconv 

# configure R.
installed.packages() |> rownames()
options('repos' = c(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN'))
install.packages('IRkernel')
IRkernel::installspec(
  name = 'omicverse.r.torch-260',
  displayname = 'omicverse.r43 (torch-2-6-0)',
  verbose = TRUE
)
# installed kernelspec omicverse.r in 
# /home/data/yangzhen/.local/share/jupyter/kernels/omicverse.r

# for r 4.3, you should install matrix 1.6-5.
install.packages('~/shared/r-pkg/matrix-1-6-5.tar.gz', repos = NULL, Ncpus = 100, type = 'source')
install.packages('~/shared/r-pkg/mass-7-3-60.tar.gz', repos = NULL, Ncpus = 100, type = 'source')
# bioconductor version of fgsea contain a dependency problem. use the development branch
install.packages('tidyverse', Ncpus = 100)
# tidyverse fails with unfound libiconv, however, this is already installed
# a workaround solution goes here.
withr::with_makevars(
  c(PKG_LIBS = '-liconv'), 
  install.packages(c('haven', 'readxl')), 
  assignment = '+='
)
BiocManager::install(c('devtools', 'remotes'), ask = FALSE, update = FALSE, Ncpus = 100)
install_github("ctlab/fgsea") 

install.packages('BiocManager', ask = FALSE, update = FALSE)
options(BioC_mirror = 'https://mirrors.westlake.edu.cn/bioconductor')
BiocManager::install('3.20')
required = c('Seurat', 'SingleR', 'clusterProfiler', 'crayon', 'tibble')
BiocManager::install(required, ask = FALSE, update = FALSE, Ncpus = 100)

options('repos' = c(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN'))
options(BioC_mirror = 'https://mirrors.westlake.edu.cn/bioconductor')

# installing jags, and infercnv
mamba install jags -c conda-forge
mamba install r-cairo
BiocManager::install(c('ggrastr', 'scater'), Ncpus = 100)
install_github ('immunogenomics/presto')
devtools::install_github("sqjin/CellChat")
mamba install hdf5
install.packages('hdf5r')
remotes::install_github("mojaveazure/seurat-disk")
