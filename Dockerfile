ARG APP=libstaden-read
ARG BUILDPLATFORM=linux/amd64
ARG TARGETARCH=amd64

FROM --platform=$BUILDPLATFORM debian:stable-slim as core

ARG TARGETARCH
ARG APP

RUN dpkg --add-architecture $TARGETARCH \
    && apt-get update \
    && apt-get install -y \      
      zlib1g-dev libbz2-dev liblzma-dev 
    



FROM core as builder

RUN apt-get install -y \
	gcc \
    wget git make \
    && rm -rf /var/lib/apt/lists/*

FROM builder as compile

RUN wget https://github.com/jkbonfield/io_lib/releases/download/io_lib-1-15-0/io_lib-1.15.0.tar.gz \
    && tar xzvf io_lib-1.15.0.tar.gz \
    && cd io_lib-1.15.0 \
    && ./configure \
    && make && make install

#FROM core as base
#
#RUN apt-get update && apt-get install -y \
#      gawk bash grep libstdc++6 libgomp1 libatomic1 zlib1g libbz2-1.0 wget tar less \
#    && rm -rf /var/lib/apt/lists/*
#
#
#FROM base as final
#
#ARG TARGETARCH
#ARG APP
#
#ARG final_path_entrypoint=/usr/local/bin/mmseqs
#
#COPY --from=builder /opt/build/${APP}_arch /opt/build/${APP}_sse2 /opt/build/${APP}_sse41 /opt/build/${APP}_avx2 /usr/local/bin/
#ADD util/${APP}_wrapper.sh ${final_path_entrypoint}
#RUN if [ "$TARGETARCH" = "arm64" ]; then \
#	rm -f ${final_path_entrypoint}; \
#	ln -s /usr/local/bin/${APP}_arch ${final_path_entrypoint}; \
#	fi
#
#ENTRYPOINT [ "/usr/local/bin/mmseqs" ]
