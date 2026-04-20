# 🐳 Guia de comandos Docker

Guia rápido com comandos essenciais do Docker para uso geral, incluindo gerenciamento de containers, imagens, volumes e uma seção dedicada ao MongoDB. 

Ao final do guia há uma tabela com a diferença entre Docker e máquinas virtuais.

---

## Imagens

- Uma imagem Docker é um template imutável que contém tudo necessário para rodar uma aplicação: código, dependências, bibliotecas e sistema base.
- Ela funciona como um “snapshot” pronto para execução.
- Imagens são construídas a partir de um Dockerfile.
- Podem ser versionadas (tags) e compartilhadas via registries.
- Servem como base para criar containers.

### Listar imagens

```bash
docker images
```
Lista todas as imagens disponíveis localmente, incluindo nome, tag e tamanho.

### Remover imagem

```bash
docker rmi <image_name>
```
Remove uma imagem específica do sistema.

Forçar remoção:
```bash
docker rmi -f <image_name>
```
Remove a imagem mesmo se estiver sendo usada por algum container.

### Construir imagem

```bash
docker build -t <image_name> .
```
Cria uma imagem a partir de um Dockerfile no diretório atual.

Sem cache:

```bash
docker build --no-cache -t <image_name> .
```
Força a reconstrução completa da imagem, ignorando camadas em cache.

### Exemplo de arquivo Dockerfile

Extrutura do projeto:

```
blast_project/
├── Dockerfile
├── requirements.txt
├── analyze_blast.py
├── query.fasta
└── pdb_sequences.fasta
```

Dockerfile:

```
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /app

# Instala BLAST+, Python e pip
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Copia os arquivos do projeto
COPY requirements.txt .
COPY analyze_blast.py .
COPY query.fasta .
COPY pdb_sequences.fasta .

# Instala dependências Python
RUN pip3 install --no-cache-dir -r requirements.txt

# Cria o banco BLAST a partir do FASTA do PDB
RUN makeblastdb -in pdb_sequences.fasta -dbtype prot -out pdb_db

# Comando padrão:
# 1) roda blastp
# 2) salva saída tabular
# 3) executa a análise em Python
CMD blastp \
    -query query.fasta \
    -db pdb_db \
    -out blast_results.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    && python3 analyze_blast.py
```

Como construir a imagem:

```
cd blast_project/
```

```bash
docker build -t blast-example .
```

---

## Containers

- Um container é uma instância em execução de uma imagem.
- Ele é isolado, leve e compartilha o kernel do sistema operacional.
- Permite rodar aplicações de forma consistente em qualquer ambiente.
- Pode ser iniciado, parado, removido e replicado facilmente.
- É efêmero por padrão (dados podem ser perdidos sem volume).

### Listar containers em execução

```bash
docker ps
```
Mostra apenas os containers ativos no momento.

### Listar todos (incluindo parados)

```bash
docker ps -a
```
Exibe todos os containers, independentemente do estado.

### Executar container interativo

```bash
docker run -it --name <container_name> <image_name>
```
Inicia um container com terminal interativo para execução manual.

### Executar em background

```bash
docker run -d --name <container_name> <image_name>
```
Inicia um container em segundo plano (modo daemon).

### Acessar container em execução

```bash
docker exec -it <container_name> bash
```
Abre um terminal dentro de um container já em execução. Para sair do container escreva "exit".

### Parar container

```bash
docker stop <container_name>
```
Encerra a execução de um container de forma segura.


Utilizando o exemplo de Dockerfile, para executar o container, use:

```bash
docker run -it --name blast_container blast-example
```

Para acessar o container em execução:

```bash
docker exec -it blast_container bash
```

### Remover container

```bash
docker rm <container_name>
```
Remove um container que já está parado.

Forçar remoção:

```bash
docker rm -f <container_name>
```
Para e remove um container ativo.

---

### Como executar um container:

Utilizando o exemplo anterior de criar uma imagem a partir de um arquivo Dockerfile, utiliza-se então esse comando para executar um container utilizando a imagem criada.

```bash
docker run -it --name blast_container blast-example
```

---

## Volumes

- Um volume é um mecanismo de persistência de dados no Docker.
- Permite armazenar dados fora do ciclo de vida do container.
- É ideal para bancos de dados e arquivos importantes.
- Pode ser compartilhado entre múltiplos containers.
- Fica armazenado no host, gerenciado pelo Docker.

### Listar volumes

```bash
docker volume ls
```
Lista todos os volumes disponíveis para persistência de dados.

### Inspecionar volume

```bash
docker volume inspect <volume_name>
```
Mostra detalhes técnicos do volume, como localização no host.

### Remover volume

```bash
docker volume rm <volume_name>
```
Remove um volume não utilizado.

### Ver volumes de um container

```bash
docker inspect <container_name>
```
Exibe informações completas do container, incluindo volumes montados.

---

## Limpeza do sistema

### Remover containers parados

```bash
docker container prune
```
Remove todos os containers que não estão em execução.

### Limpeza geral (containers, redes, volumes órfãos)

```bash
docker system prune --volumes
```
Remove recursos não utilizados para liberar espaço.

### Limpeza completa (inclui imagens)

```bash
docker system prune -a --volumes
```
Remove tudo que não está sendo usado, incluindo imagens.
**Nota**: Use com cuidado — pode remover recursos importantes.

---

## Docker Compose

- O Docker Compose é uma ferramenta para definir e gerenciar múltiplos containers.
- Utiliza um arquivo docker-compose.yml para descrever serviços, redes e volumes.
- Permite subir toda a aplicação com um único comando (docker compose up).
- Facilita ambientes complexos (ex: API + banco de dados + frontend).
- É essencial para desenvolvimento e orquestração local.

### Subir serviços

```bash
docker compose up
```
Inicia os serviços definidos no arquivo docker-compose.yml.

### Build + execução em background

```bash
docker compose up --build -d
```
Reconstrói as imagens e inicia os serviços em segundo plano.

### Parar serviços

```bash
docker compose down
```
Para e remove os containers criados pelo Compose.


### Rebuild sem cache

```bash
docker compose build --no-cache
```
Reconstrói as imagens ignorando o cache.

### Exemplo de arquivo `docker-compose.yml`:

```
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

**service app**

  Container da aplicação principal.

- build: .     constrói imagem a partir do Dockerfile
- ports        expõe a porta 8000 no host
- volumes      permite editar código localmente
- depends_on   garante que MongoDB inicie antes
- environment  variável de conexão com banco

**service mongodb**

  Container do banco de dados.

- usa imagem oficial do MongoDB
- porta padrão 27017
- volume para persistência de dados

**volumes**

- Define armazenamento persistente gerenciado pelo Docker.
- Mesmo que o container seja removido, os dados continuam salvos.

---

# MongoDB com Docker

## Executar MongoDB

```bash
docker run -d --name mongodb -p 27017:27017 -v mongodb_data:/data/db mongo:latest
```
Inicia um container MongoDB com persistência de dados em volume.

---

## Backup

### Backup a partir do host

```bash
mongodump --uri="mongodb://localhost:27017/<database>" --out ./backup
```
Cria um backup do banco diretamente no host.

### Backup dentro do container

```bash
docker exec <container_name> mongodump --uri="mongodb://localhost:27017/<database>" --out /backup
```
Executa o backup internamente no container.

### Copiar backup para o host

```bash
docker cp <container_name>:/backup ./backup_local
```
Transfere os arquivos de backup do container para o host.

### Ajustar permissões

```bash
sudo chown -R $USER:$USER backup backup_local
```
Garante que o usuário atual tenha acesso aos arquivos de backup.

---

## Restauração

### Copiar backup para o container

```bash
docker cp ./backup <container_name>:/backup
```
Move os arquivos de backup do host para o container.

### Executar restore

```bash
docker exec -it <container_name>   mongorestore --db <database> /backup
```
Restaura o banco de dados a partir do backup.

---

## Acessar MongoDB

```bash
docker exec -it <container_name> mongosh
```
Abre o shell interativo do MongoDB dentro do container.


## Docker vs Máquina Virtual

| Aspecto             | Docker (Containers)                         | Máquina Virtual (VirtualBox, VMware)       |
| ------------------- | ---------------------------------------------- | ---------------------------------------------- |
| O que é             | Ambiente isolado leve para executar aplicações | Sistema operacional completo virtualizado      |
| Sistema operacional | Compartilha o kernel do host                   | Cada VM possui seu próprio sistema operacional |
| Peso                | Leve (MB a poucos GB)                          | Pesado (vários GB)                             |
| Performance         | Próxima do nativo                              | Overhead maior                                 |
| Isolamento          | Isolamento de processos                        | Isolamento completo de hardware                |
| Uso típico          | Deploy de apps, APIs, pipelines                | Testar sistemas operacionais completos         |
| Consumo de RAM      | Baixo                                          | Alto                                           |
