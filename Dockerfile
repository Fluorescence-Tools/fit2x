# CLion remote docker environment (How to build docker container, run and stop it)
#
# Build and run:
#
#   docker build -t fit2x:0.0.5 -f Dockerfile .
#   docker run -d --cap-add sys_ptrace -p127.0.0.1:2222:22 --name fit2x_remote_env fit2x:0.0.5
#   ssh-keygen -f "$HOME/.ssh/known_hosts" -R "127.0.0.1:2222"
#
# stop:
#   docker stop fit2x_remote_env
#
# ssh credentials (test user):
#   user@password
FROM continuumio/miniconda3

SHELL ["/bin/bash", "-c"]
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update \
  && apt-get install -y ssh \
      build-essential \
      gcc \
      g++ \
      gdb \
      clang \
      rsync \
      gdb \
      tar \
      less \
      wget \
  && apt-get clean

RUN conda update -y --all

# Update the base conda environment:
COPY .condarc /root/
COPY environment.yml /root/
RUN conda env update --name base --file /root/environment.yml

# Manage ssh/sftp access
RUN ( \
    echo 'LogLevel DEBUG2'; \
    echo 'PermitRootLogin yes'; \
    echo 'PasswordAuthentication yes'; \
    echo 'Subsystem sftp /usr/lib/openssh/sftp-server'; \
  ) > /etc/ssh/sshd_config_test_clion \
  && mkdir /run/sshd

RUN useradd -m user \
  && yes password | passwd user

# create an executable cmake that activates the conda environemnt, as clion does not
# currently support env variables.
#RUN mkdir /opt/bin && echo -e '#!/usr/bin/env bash\nsource /opt/conda/etc/profile.d/conda.sh\nconda activate chinet\ncmake $@' >> /opt/bin/cmake && chmod +x /opt/bin/cmake

CMD ["/usr/sbin/sshd", "-D", "-e", "-f", "/etc/ssh/sshd_config_test_clion"]