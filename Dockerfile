ARG ALPINE_VERSION=3.12

FROM alpine:${ALPINE_VERSION} AS alpine-build

RUN apk update && apk add git cmake build-base boost-dev boost-static \
  zlib-dev zlib-static bzip2-dev bzip2-static

WORKDIR /build

COPY ./ /build/fastxio/

WORKDIR /build/fastxio/build

RUN cmake -DBUILD_APPS=ON .. && make

FROM alpine:${ALPINE_VERSION}

RUN apk add --no-cache bash

RUN apk add --no-cache zlib libbz2 libstdc++ libgcc boost-program_options libgomp

WORKDIR / 
COPY --from=alpine-build /build/fastxio/build/apps/ta_sites/TASites /usr/local/bin/
COPY --from=alpine-build /build/fastxio/build/apps/contaminant_filter/ContaminantFilter /usr/local/bin/
COPY --from=alpine-build /build/fastxio/build/apps/extr_prom/ExtractRegRegions /usr/local/bin/
COPY --from=alpine-build /build/fastxio/build/apps/kmer_map/KmerMap /usr/local/bin/
COPY --from=alpine-build /build/fastxio/build/libfastxio.so /usr/local/lib
COPY --from=alpine-build /build/fastxio/build/external/libbs/libbs.so /usr/local/lib
ENTRYPOINT ["/usr/local/bin/KmerMap"]

LABEL maintainer='github.com/bschiffthaler'
LABEL software.version=devel
