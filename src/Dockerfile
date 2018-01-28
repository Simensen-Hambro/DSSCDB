FROM continuumio/anaconda3

RUN apt-get update
RUN apt-get install build-essential zip postgresql-9.4 -y


RUN conda install -c openbabel openbabel -y
RUN conda install -c rdkit rdkit -y

RUN mkdir /code/
COPY requirements.txt /code/requirements.txt
COPY requirements_prod.txt /code/requirements_prod.txt

RUN pip install -r /code/requirements.txt
RUN pip install -r /code/requirements_prod.txt

ENTRYPOINT [ "/usr/bin/tini", "--" ]
COPY entrypoint.sh /

RUN chmod +x /entrypoint.sh
CMD ["bash", "entrypoint.sh"]
