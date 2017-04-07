FROM martsime/dsscdb:latest

RUN apt-get update
RUN apt-get install build-essential -y

ENV PATH /opt/conda/envs/rdkitenv/bin:$PATH

RUN mkdir /code/
COPY requirements.txt /code/requirements.txt
COPY requirements_prod.txt /code/requirements_prod.txt

RUN pip install -r /code/requirements.txt
RUN pip install -r /code/requirements_prod.txt

COPY entrypoint.sh /
RUN chmod +x /entrypoint.sh
