FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends g++ gcc make python3 python3-pip python3-venv curl ca-certificates bzip2 gzip git util-linux && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /venv && /venv/bin/pip install --no-cache-dir numpy scipy==1.15.3 highspy==1.15.1
ENV PATH=/venv/bin:$PATH
WORKDIR /app
COPY . /app
CMD ["bash", "/app/run_all.sh"]
