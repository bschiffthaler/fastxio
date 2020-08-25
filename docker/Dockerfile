FROM alpine:3.12 as alpine-build

RUN apk update && apk add git cmake build-base boost-dev boost-static \
  zlib-dev zlib-static bzip2-dev bzip2-static

WORKDIR /build

RUN git clone https://github.com/bschiffthaler/fastxio && \
    git checkout devel

WORKDIR /build/fastxio/build

RUN cmake -DBUILD_APPS=ON -DFASTX_STATIC_APPS=ON ..