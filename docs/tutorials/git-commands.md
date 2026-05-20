# 50 Essential Git Commands

## 1. Initial Configuration

| Command | Purpose |
| --- | --- |
| `git init` | Initialize a new Git repository in the current directory. |
| `git config --global user.name "Your Name"` | Set the global author name for commits. |
| `git config --global user.email "your.email@example.com"` | Set the global author email for commits. |
| `git config --list` | Show the current Git configuration. |
| `git config --global core.editor "editor"` | Set the default editor for commit messages, such as `nano` or `vim`. |

## 2. Repository Creation and Management

| Command | Purpose |
| --- | --- |
| `git clone <URL>` | Clone a remote repository locally. |
| `git remote add origin <URL>` | Connect the local repository to a remote GitHub repository. |
| `git remote -v` | List configured remote repositories. |
| `git remote remove origin` | Remove the remote named `origin`. |
| `git remote set-url origin <new-URL>` | Change the URL of the remote repository. |

## 3. File Management

| Command | Purpose |
| --- | --- |
| `git add <file>` | Add a specific file to the staging area. |
| `git add .` | Add all modified files to the staging area. |
| `git rm <file>` | Remove a file from the repository and stage the removal. |
| `git mv <old-file> <new-file>` | Rename or move a file. |
| `git restore <file>` | Discard uncommitted changes in a file. |
| `git restore --staged <file>` | Unstage a file while keeping local changes. |

## 4. Commits

| Command | Purpose |
| --- | --- |
| `git commit -m "Commit message"` | Create a commit from staged changes. |
| `git commit -a -m "Message"` | Automatically stage tracked files and commit them. |
| `git commit --amend` | Modify the latest commit message or contents. |
| `git log` | Show the commit history. |
| `git log --oneline` | Show one compact line per commit. |
| `git log --graph` | Show the commit history with a branch graph. |
| `git show <commit>` | Show details for a specific commit. |

## 5. Branches

| Command | Purpose |
| --- | --- |
| `git branch` | List local branches. |
| `git branch <name>` | Create a new branch. |
| `git checkout <name>` | Switch to the specified branch. |
| `git checkout -b <name>` | Create and switch to a new branch. |
| `git branch -d <name>` | Delete a local branch. |
| `git merge <name>` | Merge the specified branch into the current branch. |
| `git rebase <name>` | Reapply commits from one branch on top of another. |

## 6. GitHub Interaction

| Command | Purpose |
| --- | --- |
| `git push origin <branch>` | Push local commits to a remote branch. |
| `git push -u origin <branch>` | Set the upstream branch for future pushes. |
| `git push origin --delete <branch>` | Delete a remote branch. |
| `git pull origin <branch>` | Download and merge changes from the remote branch. |
| `git fetch origin` | Download remote changes without merging. |
| `git pull --rebase` | Pull changes and reapply local commits on top. |

## 7. Inspection and Comparison

| Command | Purpose |
| --- | --- |
| `git status` | Show the current repository state. |
| `git diff` | Show differences between modified files and the last commit. |
| `git diff --staged` | Show differences between staged files and the last commit. |
| `git blame <file>` | Show who last changed each line of a file. |

## 8. Undoing Changes

| Command | Purpose |
| --- | --- |
| `git reset --soft <commit>` | Move HEAD to a commit while keeping changes staged. |
| `git reset --hard <commit>` | Move HEAD to a commit and discard local changes. |
| `git revert <commit>` | Create a new commit that reverses another commit. |
| `git clean -f` | Remove untracked files from the working directory. |

Use destructive commands such as `git reset --hard` and `git clean -f` carefully. They can permanently remove local work.

## 9. Stash

| Command | Purpose |
| --- | --- |
| `git stash` | Temporarily save uncommitted changes. |
| `git stash list` | List saved stashes. |
| `git stash apply` | Apply the latest stash without removing it. |
| `git stash pop` | Apply the latest stash and remove it from the stash list. |
| `git stash drop` | Remove the latest stash. |
| `git stash clear` | Remove all stashed states. |
