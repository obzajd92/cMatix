# --- Stage 1: Build & Optimization Compile Layer ---
FROM ubuntu:24.04 AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY Ritz.c test.c Makefile ./

# Build execution binary footprints across optimization layers
RUN make compile_O0 compile_O2 compile_O3

# --- Stage 2: Minimal Immutable Runtime Staging Node ---
FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
    valgrind \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
COPY --from=builder /app/test_suite_O0 /workspace/test_suite_O0
COPY --from=builder /app/test_suite_O2 /workspace/test_suite_O2
COPY --from=builder /app/test_suite_O3 /workspace/test_suite_O3
COPY --from=builder /app/Ritz.c /workspace/Ritz.c

# Execute optimized solver logic as default entry barrier hook
ENTRYPOINT ["/workspace/test_suite_O3", "--benchmark", "-O3"]
