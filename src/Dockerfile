FROM martsime/dsscdb:latest

RUN apt-get update
RUN apt-get install build-essential -y

RUN mkdir /code
COPY requirements.txt /code/
COPY requirements_prod.txt /code/

ENV PATH /opt/conda/envs/rdkitenv/bin:$PATH

RUN pip install -r /code/requirements.txt
RUN pip install -r /code/requirements_prod.txt

COPY . /code/
WORKDIR /code
RUN chmod +x /code/entrypoint.sh
