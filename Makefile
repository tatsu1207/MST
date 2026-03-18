# ============================================================================
# MST-Pipeline — Makefile
# ============================================================================

.PHONY: help setup run db-reset clean check-env

ENV_NAME = mst
PYTHON = conda run -n $(ENV_NAME) python
UVICORN = conda run -n $(ENV_NAME) uvicorn

help: ## Show this help message
	@echo "MST-Pipeline — Available commands:"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'
	@echo ""

setup: ## Run the setup script
	chmod +x setup_ubuntu.sh
	./setup_ubuntu.sh

run: ## Start the application
	chmod +x run.sh
	./run.sh

db-reset: ## Reset the database (WARNING: deletes all data records, not files)
	@echo "⚠️  This will delete the database and recreate it empty."
	@read -p "Are you sure? (y/N): " confirm; \
	if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
		rm -f mst.db; \
		$(PYTHON) -c "from app.db.database import init_db; init_db()"; \
		echo "✓ Database reset."; \
	else \
		echo "Cancelled."; \
	fi

clean: ## Remove temporary files and __pycache__
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete 2>/dev/null || true
	@echo "✓ Cleaned temporary files."

clean-all: ## Remove ALL data (uploads, datasets, database)
	@echo "⚠️  This will delete ALL data including uploads and datasets."
	@read -p "Are you sure? (y/N): " confirm; \
	if [ "$$confirm" = "y" ] || [ "$$confirm" = "Y" ]; then \
		rm -rf data/uploads/* data/datasets/*; \
		rm -f mst.db; \
		echo "✓ All data removed."; \
	else \
		echo "Cancelled."; \
	fi

check-env: ## Verify environment and dependencies
	@echo "Checking environment..."
	@conda run -n $(ENV_NAME) python -c "import fastapi, dash, plotly, sqlalchemy, pandas, numpy; print('✓ Python packages OK')"
	@conda run -n $(ENV_NAME) Rscript -e 'library(dada2); message("✓ R packages OK")'
	@conda run -n ST python -c "from sourcetracker._sourcetracker import gibbs; print('✓ SourceTracker OK')"
	@echo "✓ Environment check passed."
