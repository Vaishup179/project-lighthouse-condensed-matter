FROM python:3.12.14-slim-bookworm@sha256:0f5b26b9518d002b6173fd61daad821fa340635ebfec5bba471013f9ca114579

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

RUN python -m pip install --no-cache-dir \
    numpy==2.2.5 \
    pandas==2.2.3 \
    scipy==1.15.3 \
    pytest==8.4.1

COPY data/ /app/data/

RUN mkdir -p /app/output \
    && chmod -R a+rX /app/data \
    && chmod a+rwx /app/output

CMD ["/bin/bash"]
