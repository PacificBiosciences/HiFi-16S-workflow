The files in this folder contains Singularity definition file (`.def`) for 
singularity container manual creation (if needed) and Dockerfile used for creation of Docker containers 
on Docker Hub used by both Singularity and Docker profile in the pipeline.

* To build the docker images and push to Docker Hub:
docker build -t kpinpb/pb-16s-nf-qiime:v0.2 .
docker tag kpinpb/pb-16s-nf-qiime:v0.2 kpinpb/pb-16s-nf-qiime:latest
docker push kpinpb/pb-16s-nf-qiime:v0.2
docker push kpinpb/pb-16s-nf-qiime:latest

* To build the singularity images (may need root permission):
singularity build pb-16s-pbtools.sif pb-16s-pbtools.def
