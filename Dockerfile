# ald-conformality: the chamber-to-feature ALD model with OpenFOAM and the Python feature solver in one image.
#
#   docker build -t ald-conformality . && docker run --rm ald-conformality            # the test suite
#   docker run --rm ald-conformality --quick                                          # steady flow + regression
#                                                                                     # case + Python analyses (10-20 min)
#   docker run --rm ald-conformality --full                                           # every chamber study (hours)
#   docker run --rm -it --entrypoint bash ald-conformality                            # a shell in the image
#
# Base: OpenCFD's OpenFOAM image. Tag 2512 is the newest release published on Docker Hub (there is no 2606 image;
# see paper/wp/wp7_ci_docker.md). The cases were written for v2606; nothing in them depends on a post-2512 feature.
FROM opencfd/openfoam-default:2512

ARG DEBIAN_FRONTEND=noninteractive
# The image is Ubuntu 24.04 without Python; python3 there is 3.12, the version the project is developed with.
# cmake and g++ build lib/aldkin (the kinetics library the boundary conditions link).
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 python3-venv python3-pip cmake g++ \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /work
# Python environment first (own layer, cached while the sources change)
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN python3 -m venv .venv \
 && .venv/bin/pip install --no-cache-dir --upgrade pip \
 && .venv/bin/pip install --no-cache-dir -e '.[test,crosscheck]'

# the rest of the repository (.dockerignore keeps the OpenFOAM results and the local venv out)
COPY . .

# Build the kinetics library lib/aldkin (cmake, ctest) and compile the two boundary-condition adapters against it
# into a user library inside the image (no home directory needed)
ENV ALDCONF_FOAM_USER_DIR=/work/.openfoam-user \
    ALDCONF_PYTHON=/work/.venv/bin/python \
    MPLBACKEND=Agg
RUN bash ./reproduce.sh --build

ENTRYPOINT ["bash", "/work/reproduce.sh"]
CMD ["--tests"]
