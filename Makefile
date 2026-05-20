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
LOCAL_SCRIPTS_DIR = $(shell pwd)/tools_scripts
CONTAINER_SCRIPTS_DIR = /opt/tools_scripts

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
