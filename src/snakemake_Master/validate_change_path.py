import os
import yaml

ROOT_DIR = "subworkflows"
TARGET_KEY = "output_dir"

def find_config_files():
    config_files = []
    for root, dirs, files in os.walk(ROOT_DIR):
        for file in files:
            if file == "config.yaml":
                config_files.append(os.path.join(root, file))
    return config_files

def load_yaml(file_path):
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)

def write_yaml(data, file_path):
    with open(file_path, 'w') as f:
        yaml.dump(data, f, default_flow_style=False)

def main():
    config_files = find_config_files()
    print(f"Found {len(config_files)} config.yaml files.\n")

    output_dirs = {}
    for file in config_files:
        data = load_yaml(file)
        value = data.get(TARGET_KEY, None)
        output_dirs[file] = value

    unique_values = set(output_dirs.values())

    for file, val in output_dirs.items():
        print(f"{file}: {val}")

    if len(unique_values) == 1:
        print("\n✅ All output_dir values are the same.")
    else:
        print("\n⚠️ Inconsistent output_dir values detected.")

    user_input = input("\nDo you want to set all output_dir to a unified value? (y/n): ").strip().lower()
    if user_input == 'y':
        new_value = input("Enter new unified output_dir value: ").strip()
        for file in config_files:
            data = load_yaml(file)
            data[TARGET_KEY] = new_value
            write_yaml(data, file)
        print("✅ All config.yaml files updated.")
    else:
        print("No changes made.")

if __name__ == "__main__":
    main()