# 1. Load the configuration file if it exists
ifneq ("$(wildcard .env)","")
    include .env
    export # This automatically exports all loaded variables to shell subprocesses
endif

# Variables
IMAGE_NAME = coregulation-tools
IMAGE_TAG = latest
DOCKERFILE_PATH = containers/tools.DockerFile
APPTAINER_DEF_PATH = containers/tools.def
APPTAINER_SIF_OUT = $(CONTAINER_PATH)/coregulation-tools.sif
BUILD_CONTEXT = .

LOCAL_SCRIPTS_DIR = $(shell pwd)/tools

.PHONY: install-docker-containers

docker-containers:
	@echo "Building Docker container from $(DOCKERFILE_PATH)..."
	docker build \
		-f $(DOCKERFILE_PATH) \
		-t $(IMAGE_NAME):$(IMAGE_TAG) \
		$(BUILD_CONTEXT)
	@echo "Container $(IMAGE_NAME):$(IMAGE_TAG) built successfully!"

apptainer-containers:
	@echo "Building Apptainer image..."
	# Creating the directory in case CONTAINER_PATH doesn't exist yet
	@mkdir -p $(CONTAINER_PATH)
	apptainer build $(APPTAINER_SIF_OUT) $(APPTAINER_DEF_PATH)


run-coregnet-docker:
	@echo "Running container via Docker"
	docker run --rm \
		--env-file .env \
		--platform linux/amd64 \
		-v "$(EXP_TEMP_PATH):/app/temp" \
		-v "$(DATA_PATH):/app/dataset" \
		-v "$(EXP_OUTPUT_PATH):/app/output" \
		-v "$(EXP_INPUT_PATH):/app/input" \
		-v "$(ANALYSIS_OUTPUT_PATH):/app/analysis" \
		-v "$(LOCAL_SCRIPTS_DIR):/app/script" \
		-e CODE_PATH="/app/script" \
		$(IMAGE_NAME):$(IMAGE_TAG) \
		python3 /app/script/run_coregnet.py --input $(input)

run-coregnet-apptainer:
	@echo "Running container via Apptainer"
	apptainer exec \
		--env-file .env \
		--bind "$(EXP_TEMP_PATH):/app/temp" \
		--bind "$(DATA_PATH):/app/dataset" \
		--bind "$(EXP_OUTPUT_PATH):/app/output" \
		--bind "$(EXP_INPUT_PATH):/app/input" \
		--bind "$(ANALYSIS_OUTPUT_PATH):/app/analysis" \
		--bind "$(LOCAL_SCRIPTS_DIR):/app/script" \
		--env CODE_PATH="/app/script" \
		$(APPTAINER_SIF_OUT) \
		python3 /app/script/run_coregnet.py --input $(input)