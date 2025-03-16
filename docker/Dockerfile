# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Use Ubuntu 22.04 w/CUDA 12.8 as the base image
#FROM nvidia/cuda:12.8.0-base-ubuntu22.04

RUN apt-get clean && \
    apt-get update && \
    apt-get install -y \
                       wget \
                       cmake \
                       vim \
                       g++ \
                       gdb \
                       valgrind \
                       mpich \
                       time \
                       python3\
                       pip\
                       git && \
    pip install numpy pillow 

# Build and install UPC
#mkdir upc && cd upc
#wget https://upc.lbl.gov/download/release/berkeley_upc-2022.10.0.tar.gz

#tar -xzvf berkeley_upc-2022.10.0.tar.gz
#mkdir build && cd build
#../berkeley_upc-2022.10.0/configure
#gmake ; gmake install
#RUN mkdir upcxx && cd upcxx && git clone https://bitbucket.org/berkeleylab/upcxx.git

# Copy the entrypoint file into the Docker image
COPY entrypoint.sh /entrypoint.sh

# Make the entrypoint script executable
RUN chmod +x /entrypoint.sh

# Define the entrypoint script that should be run
ENTRYPOINT ["/entrypoint.sh"]
