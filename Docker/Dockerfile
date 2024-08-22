FROM continuumio/miniconda3:main

RUN
RUN set -ex \
    && apt-get update > /dev/null \
    && apt-get install --no-install-recommends --yes \
        tini \
        > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* 

WORKDIR /
# RUN rm -rf /usr/tools/tmp
ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]