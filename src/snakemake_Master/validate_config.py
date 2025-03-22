import yaml
import jsonschema
from jsonschema import validate, ValidationError

# Load your configuration file
with open("config.yaml") as f:
    config = yaml.safe_load(f)

# Load your JSON schema file
with open("./subworkflows/snakemake_decoygenerate/workflow/config.schema.yaml") as f:
    schema = yaml.safe_load(f)

try:
    # Validate the configuration against the schema
    validate(instance=config, schema=schema)
    print("Configuration is valid!")
except ValidationError as e:
    # Print a detailed error message
    print("Validation error:", e.message)
    # The absolute_path gives a list of keys leading to the error
    print("Error path:", list(e.absolute_path))
