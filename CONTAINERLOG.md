# containerlog

Build container with `docker build -t user/reponame:tag .` given the Dockerfile is in the current directory. Then, given the repository was created at DockerHub, push it to the Hub via `docker push user/reponame:tag`. This can either be done from local machine or directly via GitPod. In the future we might add a GitHub Actions to automate that upon pushes to the `environment.yml`.

## v1.6.0
- update to Bioconductor 3.16

## v1.5.0
- updated salmon to 1.9.0 and build with current versions from the environment.yml
- update mambaforge image to 4.14.0-0

## v1.4.0
- omitted micromamba as previously used Dockerfiles now magically fail to build with non-helpful error messages -- now mambaforge again
- added [salmon binary](https://github.com/COMBINE-lab/salmon/releases/tag/v1.6.0) to PATH as the conda install for it has some issues related to dependencies with [incorrect versions](https://twitter.com/dpryan79/status/1368116490801717251?s=19), and all other software via conda
## v1.3.0
- switched to micromamba base image

## v1.2.0
- added tzdata to apt-get install to set a timezone which is not part of the stripped-down mambaforge container

## v1.1.0
- added multiqc

## v1.0.0
- first commit

