FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends g++ git ca-certificates && rm -rf /var/lib/apt/lists/*
WORKDIR /taral
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt && pip freeze | grep -i -E "numpy|scipy|highspy" > pinned_versions.txt
COPY . .
CMD ["bash", "reproduce.sh"]
