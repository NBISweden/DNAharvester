FROM python:3.11-slim-bookworm
MAINTAINER "Verena Kutschera" verena.kutschera@scilifelab.se

WORKDIR /usr/src/amber

# Install git
RUN apt-get -y update
RUN apt-get -y install git

# Clone the git repository with AMBER and install dependencies
RUN git clone https://github.com/tvandervalk/AMBER.git /usr/src/amber
RUN pip install --no-cache-dir --upgrade pysam && \
    pip install --no-cache-dir --upgrade matplotlib

COPY . .

CMD [ "python", "./AMBER" ]