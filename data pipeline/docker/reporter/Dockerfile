# syntax=docker/dockerfile:1

FROM python:3.12-slim

# ---- Configurable Quarto version ----
ARG QUARTO_VERSION=1.5.57

# If you want PDF rendering, set this to "true" at build time:
#   docker build --build-arg INSTALL_TEX=true -t myquarto .
ARG INSTALL_TEX=false

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# ---- System deps ----
# - curl/ca-certificates: download quarto
# - git: often needed if you use quarto extensions/themes
# - make: sometimes needed by project tooling
RUN apt-get update && apt-get install -y --no-install-recommends \
      curl \
      ca-certificates \
      git \
      make \
      locales \
    && rm -rf /var/lib/apt/lists/*

# ---- Install Quarto ----
# Quarto bundles Pandoc, so you usually don't need to install pandoc separately.
RUN curl -fsSL -o /tmp/quarto.deb \
      "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.deb" \
    && apt-get update && apt-get install -y --no-install-recommends /tmp/quarto.deb \
    && rm -f /tmp/quarto.deb \
    && rm -rf /var/lib/apt/lists/*

# ---- Optional: LaTeX for PDF output ----
# This is intentionally minimal-ish. If your PDF builds fail due to missing TeX packages,
# you may need a larger texlive set (texlive-full is huge).
RUN if [ "$INSTALL_TEX" = "true" ]; then \
      apt-get update && apt-get install -y --no-install-recommends \
        texlive-latex-base \
        texlive-latex-recommended \
        texlive-latex-extra \
        texlive-fonts-recommended \
        texlive-xetex \
        latexmk \
      && rm -rf /var/lib/apt/lists/* ; \
    fi

# ---- (Optional) Install Python deps ----
# Copy your requirements separately for better Docker caching
WORKDIR /app
COPY requirements.txt /app/requirements.txt
RUN if [ -f requirements.txt ]; then pip install --no-cache-dir -r requirements.txt; fi

# Copy the rest of your project (including .qmd, data, etc.)
COPY . /app

# Sanity check: show Quarto version during build logs
RUN quarto --version && python --version


COPY reporting_module/requirements.txt .
RUN pip install -r requirements.txt

# Default command: render a Quarto project or a single file
# Override with `docker run ... quarto render <something>`
CMD ["quarto", "render"]