import pipeline
from pathlib import Path

if __name__ == "__main__":
    pipeline.run_halper_pipeline(Path("config.yaml"))