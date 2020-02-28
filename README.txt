# NGS snp calling pipeline

A bash pipeline for creating a *de novo* assembly of Wild-type which is compared to mutant sequences to identifyputative SNPs and indels

## Installations 

Download the zip file and unzip in chosen directory.

```bash
unzip NGSdownstreamPipeline.zip
```

# Install dependencies

Use [Linuxbrew](https://docs.brew.sh/Homebrew-on-Linux) for installing most dependencies

```bash
brew install brewsci/bio/prokka
brew install brewsci/bio/shovill
shovill --check

brew install vt
brew install bwa
brew install snpEff
brew install SnpSift

easy_install xlsx2csv ~OR~ pip xlsx2csv

sudo apt-get update -y
sudo apt-get install -y tabix
sudo apt-get install -y bgzip
```

## Usage

```bash
bash runNGSdownstreamPipeline.sh -r [PATH_TO_R1_AND_R2_READS] -m [META_FILE.xlsx]
```


