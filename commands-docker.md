# Docker commands

## Comandos básicos

### Listar as imagens

```bash
docker images
```

### Deletar uma imagem

```bash
docker rmi name_image
```

ou para forçar a deleção

```bash
docker rmi -f name_image
```

### Listar os contêiners

```bash
docker ps
```

### Construir uma imagem (sendo o nome app do diretório de trabalho)

```bash
docker build -t app .
```

### Executar conteiner

```bash
docker run -it --name nome-do-conteiner app
```

### Acessar o shell de um contêiner Docker em execução

```bash
docker exec -it -u root id-do-conteiner bash
```

### Inicia um contêiner interativo a partir da imagem app e abre um shell (sh) dentro do contêiner

```bash
sudo docker run -it app bash
```
### Construir a imagem antes de executar o conteiner em segundo plano (flag -d)

```bash
docker-compose up --build -d
```

### Construir imagens sem usar o cache

```bash
docker compose build
```

### Listar Volumes
```bash
docker volume ls
```

### Inspecionar volume
```bash
docker volume inspect nome-do-conteiner
```

### Ver volumes montados em um container específico

```bash
docker inspect nome-do-conteiner | grep -A 10 -B 10 Mount
```

### Apagar volume
```bash
docker volume rm nome-do-conteiner
```

## Backup (MongoDB)
Backup Usando ```mongodump``` diretamente do host
```bash
sudo mongodump --uri="mongodb://172.17.0.2:27017/gh32" --out ./backup
```
O backup será salvo na pasta backup do host

Ou executar o backup dentro do próprio container MongoDB
```bash
docker exec mongodb mongodump --uri="mongodb://172.17.0.2:27017/gh32" --out /backup
```

### Copiar para o seu host
```bash
sudo docker cp mongodb:/backup ./backup_local
```
Será criado um backup na pasta ./backup

### 
Alterando a propriedade das pastas de root para o usuário
```bash
sudo chown -R $USER:$USER backup backup_local
```

### Restauração do banco

Transferindo a pasta contendo o backup para o conteiner
```bash
docker cp /caminho/do/backup mongodb:/backup
```

Acessar o conteiner
```bash
docker exec -it mongodb sh
```

Executando o mongorestore para restaurar o backup
```bash
mongorestore --db gh32 /backup
```

## Limpeza de Containers, Volumes e Redes no Docker

Comandos úteis para remover containers, volumes e redes **não utilizados**, ajudando a manter o ambiente Docker limpo e eficiente.

### 1. Remover apenas containers parados

Remove **somente** containers que não estão em execução.

```bash
docker container prune
```

Para executar sem confirmação interativa:

```bash
docker container prune -f
```

### 2. Remover todos os containers (inclusive ativos)

```bash
docker rm -f $(docker ps -aq)
```

**Atenção**: esse comando para e remove tudo.
Use apenas se quiser limpar completamente os containers.

### 3. Limpar containers, redes e volumes não utilizados

Remove:

Containers **parados**

**Redes** que não estão em uso

**Volumes órfãos** (não associados a nenhum container)

**Caches** de build

```bash
docker system prune --volumes
```

Para executar sem confirmação

```bash
docker system prune --volumes -f
```

**Não remove imagens.**
É o comando ideal para uma limpeza segura e completa de containers e recursos inativos.


### 4. Limpeza total (containers + imagens + volumes + redes)

**Remove tudo que não está sendo usado** — inclusive **imagens** que não estão associadas a nenhum container.

```bash
docker system prune -a --volumes
```
---

### Run DLKcat

Para obter a imagem

```bash
docker pull aetafur/dlkcat
```

Execute no diretório atual de trabalho
```bash
docker run -d --rm --gpus all --name dlkcat-container -v "$(pwd)":/workspace -w /workspace aetafur/dlkcat:latest
```

Para entrar no shell do container em execução:
```bash
docker exec -it dlkcat-container bash
```

Para parar o conteiner:
```bash
docker stop dlkcat-container
```
Como foi usado `--rm`, ele será removido automaticamente ao parar.

---
## Docker compose commands
> https://docs.docker.com/compose/reference/
