## Introduction
### This project is for rebuiding gene.iobio database based on [gene.model.db](https://github.com/iobio/gene.model.db/tree/master)

### Since the gff files are massive, it will take a while to process and read into the database. That said, there is space to optimize. For example, use Dask to process data parallelly.

____

## Instruction
### 1. Install 
```git clone project ```
### 2. Create your virtualenv
``` python3 -m venv myenv ```
### 3. Install packages and dependencies
``` pip3 install -r requirements.txt```
### 4. Run
``` Python3 data_downloader.py``` <br>
``` Python3 database_builder.py```

