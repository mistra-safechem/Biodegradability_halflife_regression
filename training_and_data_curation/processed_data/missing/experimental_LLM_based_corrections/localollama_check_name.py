"""Script to correct missing chemical names using local LLM model.
This works in principal but returns too many incorrect results to be useful and thus would require manual checking anyway.
Nice proof of concept though.
Alexander Minidis, 2025-10-23

Remarks in proof: large LLMs are much better at understanding and correcting chemical names than smaller ones (as experienced in different projects).
Requires though high memory local system or cloud llm access - question of token pricing (and optimizing prompts which they aren't at this stage).

"""

# Check if ollama is installed, if not, install with 'uv add ollama'
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

if importlib.util.find_spec("ollama") is None:
    print("ollama not found, installing with 'uv add ollama'...")
    subprocess.check_call([sys.executable, "-m", "uv", "add", "ollama"])

import ollama
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def main():
    script_dir = Path(__file__).parent.resolve()
    missing_names = script_dir.parent / "missing_names_combo_all_compartments.csv"
    found_names = script_dir / "missing_names_combo_all_compartments_llm_test.csv"

    missing_smiles = []
    for line in open(missing_names, "r"):
        line = line.strip()
        if line:
            missing_smiles.append(line)

    # two different modeles and prompts, choose one
    model = "gemma4:e4b"  # works to some extent, but not very good with chemistry names
    prompt_template = """You are a chemical expert proficient in analyzing generic substance names to find faults and correct them. 
    Here is a generic substance name: {input_text}. Please correct it if necessary. Then provide the IUPAC name of the substance.
    Output format: JSON '{{"input": "{input_text}", "output": "IUPAC name"}}', with the 'IUPAC name' being the single result, nothing else.
    Start your response directly with the processed text and no acknowledgements about my questions and no explanations required.
    """
    # model_kale_LM = "hf.co/mradermacher/Llama3-KALE-LM-Chem-1.5-8B-GGUF:Q8_0"  # chemistry centric llm. works a bit better?
    # prompt_kale_LM_template = """correct the following generic substance name {input_text} and return input and identified IUPAC name as JSON format: '{{"input": "{input_text}", "output": "IUPAC name"}}'."""
    # model = model_kale_LM
    # prompt = prompt_kale_LM_template

    new_names = []
    print(f"Total names to correct: {len(missing_smiles)}")
    for name in missing_smiles:  # slice if you want to test with a subset
        input_text = name
        prompt = prompt_template.format(input_text=input_text)
        # prompt_kale_LM = prompt_kale_LM_template.format(input_text=input_text)

        # select model and prompt here
        response = ollama.generate(model=model, prompt=prompt, options={"temperature": 0})

        actual_response = response["response"].strip()
        # print(f"Raw response for input '{name}': {actual_response}")

        new_name = "NA"
        try:
            actual_response_json = json.loads(actual_response)
            new_name = actual_response_json.get("output", "")
        except json.JSONDecodeError:
            print(f"Warning: Invalid JSON for input '{name}': {actual_response}")
            if " is " in actual_response:
                new_name = actual_response.split(" is ")[-1].strip()
            if "IUPAC" in actual_response:
                new_name = "NA"
        print(f"input name: {name}. New name: {new_name}")
        new_names.append(new_name)

    combine_names = list(zip(missing_smiles, new_names))
    # note: the kalelm output might contain multiple names with ; which sort of messes up df output; not a big deal but still
    df = pd.DataFrame(combine_names, columns=["original_name", "corrected_name"])
    df.to_csv(found_names, index=False, header=True, sep="\t")


if __name__ == "__main__":
    main()
