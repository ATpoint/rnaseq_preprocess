# containerlog

Build container with `docker build -t user/reponame:tag .` given the Dockerfile is in the current directory. Then, given the repository was created at DockerHub, push it to the Hub via `docker push user/reponame:tag`. This can either be done from local machine or directly via GitPod. In the future we might add a GitHub Actions to automate that upon pushes to the `environment.yml`.

## v1.3.0
- switched to micromamba base image

## v1.2.0
- added tzdata to apt-get install to set a timezone which is not part of the stripped-down mambaforge container

## v1.1.0
- added multiqc

## v1.0.0
- first commit

