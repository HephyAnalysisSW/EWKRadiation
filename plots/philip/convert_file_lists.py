import os

def convert_filelist(input_path, output_path):
    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        current_base = ""
        for line in fin:
            line = line.strip()
            if not line:
                continue
            if line.endswith(":"):
                current_base = line[:-1]  # Entferne den Doppelpunkt
            elif current_base:
                full_path = os.path.join(current_base, line)
                fout.write(full_path + "\n")

def convert_all_lists(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(input_dir):
        if filename.endswith(".txt"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, filename.replace(".txt", "_cleaned.txt"))
            print(f"Converting: {filename} -> {os.path.basename(output_path)}")
            convert_filelist(input_path, output_path)

if __name__ == "__main__":
    input_dir = "file_lists"           # Ordner mit den RAW-Listen
    output_dir = "cleaned_file_lists"  # Neuer Ordner mit bereinigten Listen
    convert_all_lists(input_dir, output_dir)
    print("✅ Done converting all file lists.")
