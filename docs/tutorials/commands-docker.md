# Docker Commands

This quick guide covers essential Docker commands for general use, including containers, images, volumes, Docker Compose, and a MongoDB section.

## Images

A Docker image is an immutable template that contains everything needed to run an application: code, dependencies, libraries, and the base system. Images are built from a `Dockerfile`, can be versioned with tags, and are shared through registries.

### List Images

```bash
docker images
```

### Remove an Image

```bash
docker rmi <image_name>
```

Force removal:

```bash
docker rmi -f <image_name>
```

### Build an Image

```bash
docker build -t <image_name> .
```

Build without cache:

```bash
docker build --no-cache -t <image_name> .
```

### Example Dockerfile

Project structure:

```text
blast_project/
|-- Dockerfile
|-- requirements.txt
|-- analyze_blast.py
|-- query.fasta
`-- pdb_sequences.fasta
```

`Dockerfile`:

```dockerfile
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /app

# Install BLAST+, Python, and pip
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Copy project files
COPY requirements.txt .
COPY analyze_blast.py .
COPY query.fasta .
COPY pdb_sequences.fasta .

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Create the BLAST database from the PDB FASTA file
RUN makeblastdb -in pdb_sequences.fasta -dbtype prot -out pdb_db

# Default command:
# 1. Run blastp
# 2. Save tabular output
# 3. Run the Python analysis
CMD blastp \
    -query query.fasta \
    -db pdb_db \
    -out blast_results.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    && python3 analyze_blast.py
```

Build the image:

```bash
cd blast_project/
docker build -t blast-example .
```

## Containers

A container is a running instance of an image. It is isolated, lightweight, and shares the host operating-system kernel. Containers are ephemeral by default, so persistent data should be stored in volumes.

### List Running Containers

```bash
docker ps
```

### List All Containers

```bash
docker ps -a
```

### Run an Interactive Container

```bash
docker run -it --name <container_name> <image_name>
```

### Run in the Background

```bash
docker run -d --name <container_name> <image_name>
```

### Access a Running Container

```bash
docker exec -it <container_name> bash
```

Type `exit` to leave the container shell.

Using the BLAST example:

```bash
docker run -it --name blast_container blast-example
docker exec -it blast_container bash
```

### Stop a Container

```bash
docker stop <container_name>
```

### Remove a Container

```bash
docker rm <container_name>
```

Force removal:

```bash
docker rm -f <container_name>
```

## Volumes

A volume stores persistent data outside the container lifecycle. It is useful for databases, analysis outputs, and files that must survive container recreation.

### List Volumes

```bash
docker volume ls
```

### Inspect a Volume

```bash
docker volume inspect <volume_name>
```

### Remove a Volume

```bash
docker volume rm <volume_name>
```

### Inspect Container Mounts

```bash
docker inspect <container_name>
```

## System Cleanup

### Remove Stopped Containers

```bash
docker container prune
```

### General Cleanup

```bash
docker system prune --volumes
```

### Full Cleanup

```bash
docker system prune -a --volumes
```

This removes unused containers, networks, volumes, and images. Use it carefully because it can remove resources you still need.

## Docker Compose

Docker Compose defines and manages multi-container applications through a `docker-compose.yml` file. It is useful for local stacks such as an API, database, and frontend.

### Start Services

```bash
docker compose up
```

### Build and Run in the Background

```bash
docker compose up --build -d
```

### Stop Services

```bash
docker compose down
```

### Rebuild without Cache

```bash
docker compose build --no-cache
```

### Example `docker-compose.yml`

```yaml
services:
  app:
    build: .
    container_name: app_container
    ports:
      - "8000:8000"
    volumes:
      - .:/app
    depends_on:
      - mongodb
    environment:
      - MONGO_URI=mongodb://mongodb:27017/mydatabase
    restart: unless-stopped

  mongodb:
    image: mongo:latest
    container_name: mongodb
    ports:
      - "27017:27017"
    volumes:
      - mongodb_data:/data/db
    restart: unless-stopped

volumes:
  mongodb_data:
```

The `app` service builds the main application image, exposes port `8000`, mounts the current directory into `/app`, waits for MongoDB to start, and uses `MONGO_URI` for the database connection.

The `mongodb` service uses the official MongoDB image, exposes port `27017`, and stores data in the `mongodb_data` volume.

## MongoDB with Docker

### Run MongoDB

```bash
docker run -d --name mongodb -p 27017:27017 -v mongodb_data:/data/db mongo:latest
```

### Backup from the Host

```bash
mongodump --uri="mongodb://localhost:27017/<database>" --out ./backup
```

### Backup from inside the Container

```bash
docker exec <container_name> mongodump --uri="mongodb://localhost:27017/<database>" --out /backup
```

### Copy Backup to the Host

```bash
docker cp <container_name>:/backup ./backup_local
```

### Fix Permissions

```bash
sudo chown -R $USER:$USER backup backup_local
```

### Copy Backup to the Container

```bash
docker cp ./backup <container_name>:/backup
```

### Restore Backup

```bash
docker exec -it <container_name> mongorestore --db <database> /backup
```

### Open the MongoDB Shell

```bash
docker exec -it <container_name> mongosh
```

## Docker vs Virtual Machine

| Aspect | Docker containers | Virtual machine |
| --- | --- | --- |
| Purpose | Lightweight isolated environment for applications | Full virtualized operating system |
| Operating system | Shares the host kernel | Each VM has its own operating system |
| Size | Lightweight, from MB to a few GB | Heavy, often several GB |
| Performance | Close to native | Higher overhead |
| Isolation | Process-level isolation | Full hardware-level isolation |
| Common use | Apps, APIs, pipelines, deployment | Testing full operating systems |
| RAM usage | Lower | Higher |
