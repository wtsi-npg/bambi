FROM debian:bookworm AS htslib_build

ARG HTSLIB_BRANCH

# hadolint ignore=DL3003
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get upgrade -y && \
    apt-get install -qq -y --no-install-recommends \
      autoconf \
      automake \
      autotools-dev \
      build-essential \
      ca-certificates \
      git \
      libbz2-dev \
      libcurl4-openssl-dev \
      libgd-dev \
      liblzma-dev \
      libssl-dev \
    && \
    export HTSLIB=/usr/local/src/htslib && \
    git clone --depth 1 --branch "${HTSLIB_BRANCH}" https://github.com/samtools/htslib.git "${HTSLIB}" && \
    cd "${HTSLIB}" && \
    git submodule update --init --recursive && \
    mkdir -p /tmp/usr/local && autoreconf -i && ./configure --prefix=/tmp/usr/local && make && make install

FROM debian:bookworm AS bambi_build

COPY --from=htslib_build /tmp/usr/local /usr/local
COPY . /usr/local/src/bambi

WORKDIR /usr/local/src/bambi

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get upgrade -y && \
    apt-get install -qq -y --no-install-recommends \
      autoconf \
      automake \
      autotools-dev \
      build-essential \
      bzip2 \
      git \
      libbz2-dev \
      libcurl4-openssl-dev \
      libtool \
      libgd-dev \
      liblzma-dev \
      libxml2-dev \
      libz-dev \
      libssl-dev \
    && \
    export PREFIX=/tmp/usr/local && \
    autoreconf -fi && \
    CPPFLAGS="-I$PREFIX/include" LDFLAGS="-Wl,-rpath,$PREFIX/lib -L$PREFIX/lib" \
            ./configure --prefix="$PREFIX" --with-htslib="$PREFIX" && \
    make CPPFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib" && \
    make install prefix="$PREFIX"

FROM debian:bookworm AS bambi

COPY --from=htslib_build /tmp/usr/local /usr/local
COPY --from=bambi_build /tmp/usr/local /usr/local

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get upgrade -y && \
    apt-get install -qq -y --no-install-recommends \
      libbz2-1.0 \
      libcurl4 \
      libgd3 \
      liblzma5 \
      libxml2 \
      libz1 \
      libssl3 \
    && \
    apt-get remove -q -y unattended-upgrades && \
    apt-get autoremove -q -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]
