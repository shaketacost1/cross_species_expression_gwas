# Makefile for building and running the frontend Docker container locally

IMAGE_NAME=frontend-app
CONTAINER_NAME=frontend-app-container
DOCKERFILE_PATH=Dockerfile
APP_DIR=.
PORT=8080
INTERNAL_PORT=80

.PHONY: build run stop clean logs status test health

build:
	docker build --platform linux/amd64 -t $(IMAGE_NAME) -f $(DOCKERFILE_PATH) $(APP_DIR)

run:
	docker run --rm -d --name $(CONTAINER_NAME) -p $(PORT):$(INTERNAL_PORT) $(IMAGE_NAME)

stop:
	docker stop $(CONTAINER_NAME) || echo "Container $(CONTAINER_NAME) not running"

clean:
	docker rm -f $(CONTAINER_NAME) 2>/dev/null || true
	docker rmi -f $(IMAGE_NAME) 2>/dev/null || true

logs:
	docker logs -f $(CONTAINER_NAME)

status:
	@echo "=== Container Status ==="
	docker ps -a --filter name=$(CONTAINER_NAME)
	@echo ""
	@echo "=== Image Info ==="
	docker images $(IMAGE_NAME)

test:
	@echo "Testing connection to http://localhost:$(PORT)..."
	@curl -s -o /dev/null -w "HTTP Status: %{http_code}\n" http://localhost:$(PORT) || echo "Connection failed"

health:
	@echo "=== Health Check ==="
	docker inspect $(CONTAINER_NAME) --format='{{.State.Health.Status}}' 2>/dev/null || echo "Container not found or no health check"

restart: clean build run

debug:
	docker run --rm --name $(CONTAINER_NAME) -p $(PORT):80 $(IMAGE_NAME)