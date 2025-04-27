import yaml
from pathlib import Path

def update_config(config_path: Path, config_entries: dict[str, Path], backup_suffix: str = ".backup") -> None:
    """
    Update the config file with the given config entries.
    """
    # Create a backup of the original config file
    backup_path = Path(f"{config_path}{backup_suffix}")
    with open(config_path, "r") as src:
        with open(backup_path, "w") as dst:
            dst.write(src.read())

    # Read the config file
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Update the config file with the given config entries
    existing_keys = []
    for key, value in config_entries.items():
        if key not in config:
            config[key] = value.as_posix()
        else:
            existing_keys.append(key)

    if config["warn_on_existing_keys"]:
        print(f"The following keys have been specified in the config file: {existing_keys}")
    
    # Write the updated config file
    with open(config_path, "w") as f:
        yaml.dump(config, f)

