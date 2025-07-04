# 50 comandos essenciais do Git

### 1. Configuração Inicial

git init - Inicializa um novo repositório Git em um diretório.

git config --global user.name "Seu Nome" - Define o nome do usuário para commits globalmente.

git config --global user.email "seu.email@exemplo.com" - Define o e-mail do usuário para commits globalmente.

git config --list - Exibe as configurações do Git.

git config --global core.editor "editor" - Define o editor padrão para mensagens de commit (ex.: "nano" ou "vim").


### 2. Criação e Gerenciamento de Repositórios
git clone <URL> - Clona um repositório remoto para o diretório local.

git remote add origin <URL> - Conecta o repositório local a um repositório remoto no GitHub.

git remote -v - Lista os repositórios remotos configurados.

git remote remove origin - Remove a conexão com o repositório remoto.

git remote set-url origin <nova-URL> - Altera a URL do repositório remoto.


### 3. Gerenciamento de Arquivos
git add <arquivo> - Adiciona um arquivo específico ao stage (área de preparação).

git add . - Adiciona todos os arquivos modificados ao stage.

git rm <arquivo> - Remove um arquivo do repositório e do stage.

git mv <arquivo-antigo> <arquivo-novo> - Renomeia ou move um arquivo.

git restore <arquivo> - Descarta alterações em um arquivo não commitado.

git restore --staged <arquivo> - Remove um arquivo do stage, mantendo as alterações locais.


### 4. Commits
git commit -m "Mensagem do commit" - Cria um commit com as alterações no stage.

git commit -a -m "Mensagem" - Adiciona e faz commit de arquivos rastreados automaticamente.

git commit --amend - Modifica o último commit (mensagem ou arquivos).

git log - Exibe o histórico de commits.

git log --oneline - Exibe o histórico de commits em uma linha por commit.

git log --graph - Mostra o histórico de commits com um gráfico de branches.

git show <commit> - Exibe detalhes de um commit específico.


### 5. Branches
git branch - Lista todas as branches locais.

git branch <nome> - Cria uma nova branch.

git checkout <nome> - Muda para a branch especificada.

git checkout -b <nome> - Cria e muda para uma nova branch.

git branch -d <nome> - Deleta uma branch local.

git merge <nome> - Mescla a branch especificada na branch atual.

git rebase <nome> - Reaplica commits de uma branch sobre outra.


### 6. Interação com o GitHub
git push origin <branch> - Envia commits locais para a branch remota no GitHub.

git push -u origin <branch> - Define a branch remota como padrão para pushes futuros.

git push origin --delete <branch> - Deleta uma branch remota.

git pull origin <branch> - Baixa e mescla alterações do repositório remoto.

git fetch origin - Baixa alterações do repositório remoto sem mesclar.

git pull --rebase - Faz pull e reaplica commits locais sobre os remotos.


### 7. Inspeção e Comparação
git status - Mostra o estado atual do repositório (arquivos modificados, staged, etc.).

git diff - Exibe as diferenças entre arquivos modificados e o último commit.

git diff --staged - Mostra as diferenças entre o stage e o último commit.

git blame <arquivo> - Mostra quem modificou cada linha de um arquivo.


### 8. Desfazer Alterações
git reset --soft <commit> - Desfaz um commit, mantendo alterações no stage.

git reset --hard <commit> - Desfaz um commit e descarta todas as alterações.

git revert <commit> - Cria um novo commit que desfaz as alterações de um commit específico.

git clean -f - Remove arquivos não rastreados do diretório de trabalho.


### 9. Stash (Armazenamento Temporário)
git stash - Salva alterações não commitadas temporariamente.

git stash list - Lista todos os stashes salvos.

git stash apply - Aplica o último stash sem removê-lo.

git stash pop - Aplica o último stash e o remove.

git stash drop - Remove o último stash.

git stash clear - Remove all the stashed states.
