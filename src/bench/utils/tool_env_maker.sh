#!/bin/bash


# Function to create environment and install tool
create_env() {
    local env_name=$1
    local package=$2
    local channel=${3:-bioconda}  # default to bioconda if no channel specified

    echo "Creating environment: $env_name with package: $package"
    
    # Remove environment if it exists (ignore errors if it doesn't exist)
    micromamba env remove -n $env_name 2>/dev/null || true
    
    # Create new environment with only the specific tool
    micromamba create -n $env_name -c $channel $package -y
}


tools=("bowtie" "bowtie2" "minimap2" "bbmap" "strobealign" "blast" "mmseqs2" "mummer4" "lexicmap" "bwa" "sassy" "x-mapper")

# Create individual environments for each tool
for tool in "${tools[@]}"; do
    if [ "$tool" != "sassy" ]; then
        create_env "${tool}_env" "$tool"
    fi
done

# Create sassy environment separately (empty environment for Rust tool)
echo "Creating sassy environment..."
micromamba env remove -n sassy_env 2>/dev/null || true
micromamba create -n sassy_env -y

# For spacer-containment, assuming it's available through a specific channel
# Adjust the channel and package name as needed
create_env "spacer_containment_env" "spacer-containment" 

echo "All environments created successfully!"




# ###### 
# # mmseqs2 16.747c6 (current bioconda version) results in frequent crashes. as it is borken, we will install the latest version from github
# micromamba activate mmseqs2_env
# wget https://mmseqs2.com/latest/mmseqs2-linux-avx2.tar.gz 
# tar -xvzf mmseqs2-linux-avx2.tar.gz
# mv ./mmseqs2/bin/mmseqs2 $CONDA_PREFIX/bin/mmseqs2
# chmod +x $CONDA_PREFIX/bin/mmseqs2
# rm -rf mmseqs2-linux-avx2.tar.gz mmseqs2
# echo "#####\n mmseqs2 version: \n" >> tool_ tool_configs/tool_versions.txt
# mmseqs2 version >> tool_ tool_configs/tool_versions.txt

###### 
# mummer 4 from main branch has bugs in SAM format output (see https://github.com/mummer4/mummer/issues/24)
# we will install the latest version of the develop branch (https://github.com/mummer4/mummer/tree/develop)
micromamba activate mummer4_env
# install dependencies
micromamba install gcc make yaggo -y
wget https://github.com/mummer4/mummer/archive/refs/heads/develop.zip # https://github.com/mummer4/mummer/releases/download/v4.0.1/mummer-4.0.1.tar.gz testing now
rm -rf mummer-develop
unzip develop.zip
cd mummer-develop
rm $CONDA_PREFIX/bin/mummer -rf
autoreconf -fi
./configure --prefix=$CONDA_PREFIX/bin/mummer
make  -j 6
make install
cd ..
rm -rf mummer-develop develop.zip
cp $CONDA_PREFIX/bin/mummer/bin/nucmer $CONDA_PREFIX/bin/nucmer
chmod +x $CONDA_PREFIX/bin/nucmer
echo "#####\n Mummer version: \n" >> tool_ tool_configs/tool_versions.txt
nucmer --version >> tool_ tool_configs/tool_versions.txt

# ######
# # hisat2
# # create_env "hisat2_env" "hisat2"
# micromamba activate hisat2_env
# echo "#####\n Hisat2 version: \n" >> tool_ tool_configs/tool_versions.txt
# hisat2 --version >> tool_ tool_configs/tool_versions.txt
# # wget https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz
# # tar -xvzf v2.2.1.tar.gz
# # cd hisat2-2.2.1
# # ./install_hisat2.sh
# # cd ..
# # rm -rf v2.2.1.tar.gz hisat2-2.2.1

# ######
# # Bwa   
# # create_env "bwa_env" "bwa"
# micromamba activate bwa_env
# echo "#####\n Bwa version: \n" >> tool_ tool_configs/tool_versions.txt
# bwa 2>&1 |  head -n 3 |tail -n 1 >> tool_ tool_configs/tool_versions.txt





######
# Tools built from source on a "benchy" environment (see benchy_env.sh for the mamba environment details)
# the pixi default environment includes maturin so rust  and cargo already installed
pixi shell
# Sassy - Rust tool that needs to be built from source
echo "sassy/" >> .gitignore # only want to compile it
git clone https://github.com/RagnarGrootKoerkamp/sassy.git
cd sassy
cargo build --release # --features "python,scalar"
# Copy the binary to the conda environment
cp target/release/sassy $CONDA_PREFIX/bin/sassy
chmod +x $CONDA_PREFIX/bin/sassy
cd ..
rm -rf sassy




# Print versions of installed tools
echo "Checking installed versions:"
echo -e "tool_name\tversion" > tool_configs/tool_versions.tsv

# Function to check if a tool is installed and print the version
check_tool_version() {
    tool_env=$1
    version_cmd=$2
    tool_name=$(basename "$tool_env" _env)
    # Run the version command in the environment and get the first line of output
    version=$(micromamba run -n "$tool_env" bash -c "$version_cmd" 2>&1 | head -n 1 | tr -d '\r')
    # If version is empty, set to "N/A"
    if [ -z "$version" ]; then
        version="N/A"
    fi
    echo -e "${tool_name}\t${version}" >> tool_configs/tool_versions.tsv
}
# tools that write to stdout and not stderr
check_tool_version "bowtie_env" "bowtie-align-s --version | head -n 1 | tr -d '\r'"
check_tool_version "bowtie2_env" "bowtie2-align-s --version | head -n 1 | tr -d '\r'"
check_tool_version "minimap2_env" "minimap2 --version"
check_tool_version "blast_env" "blastn -version" # blast executable is called blastn
check_tool_version "strobealign_env" "strobealign --version"
check_tool_version "mmseqs2_env" "mmseqs version" # mmseqs executable is called mmseqs, without th 2
# tools that write to stderr and not stdout 
check_tool_version "bbmap_env" "bbmap.sh --version 2>&1 | head -n 1 | tr -d '\r'" # bbmap execs are wrapped by shell scripts (.sh), without the dot
check_tool_version "x-mapper_env" "x-mapper --version"
check_tool_version "benchy" "sassy --version"
