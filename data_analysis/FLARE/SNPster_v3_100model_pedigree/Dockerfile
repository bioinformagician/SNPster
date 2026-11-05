FROM snpster-ancestry-benchmark:v1

ARG PEDSIM_COMMIT=44d0c7abf6206479d5d6df926ecca086a8bbee32

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends libboost-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/williamslab/ped-sim.git /opt/ped-sim \
    && git -C /opt/ped-sim checkout --detach "${PEDSIM_COMMIT}" \
    && make -C /opt/ped-sim -j2 \
    && ln -s /opt/ped-sim/ped-sim /usr/local/bin/ped-sim

COPY requirements.txt /tmp/pedigree-requirements.txt
RUN pip3 install --no-cache-dir -r /tmp/pedigree-requirements.txt

WORKDIR /workspace
