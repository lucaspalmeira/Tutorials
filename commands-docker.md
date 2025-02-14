## Docker compose commands
> https://docs.docker.com/compose/reference/

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
docker exec -it -u root id-do-conteiner sh
```

### Inicia um contêiner interativo a partir da imagem app e abre um shell (sh) dentro do contêiner

```bash
sudo docker run -it app sh
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
docker volume inspect mongo-data
```

### Apagar volume
```bash
docker volume rm mongo-data
```

### Backup
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
