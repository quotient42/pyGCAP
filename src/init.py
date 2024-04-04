import os
import pandas as pd

project_info = {
    'Project_Name': [],
    'Root_Dir': [],
    'Input_Dir': [],
    'Output_Dir': [],
    'Data_Dir': [],
    'Seqlib_Dir': []
}

root_dir = input(">> Root Directory: ")
if not os.path.exists(root_dir):
    print(f"<< Directory not found. Check your input: {root_dir}")
    exit()
elif not os.path.isdir(root_dir):
    print(f"<< Not a directory. Check your input: {root_dir}")
    exit()

project_name = input(">> Project Name: ")
if not project_name:
    project_name = "new project"

try:
    for dir_name in ['input', 'data', 'seqlib', 'output', 'output/tsv', 'output/img', 'output/genus']:
        if not os.path.isdir(os.path.join(root_dir, dir_name)):
            os.makedirs(os.path.join(root_dir, dir_name))
except Exception as e:
    print(f"<< Directory creation failed: {e}")

project_info['Project_Name'].append(project_name)
project_info['Root_Dir'].append(root_dir)
project_info['Input_Dir'].append(os.path.join(root_dir, 'input'))
project_info['Output_Dir'].append(os.path.join(root_dir, 'output'))
project_info['Data_Dir'].append(os.path.join(root_dir, 'data'))
project_info['Seqlib_Dir'].append(os.path.join(root_dir, 'seqlib'))

df_project_info = pd.DataFrame(project_info)

df_project_info.to_csv('project_info.tsv', sep='\t', index=False)

print(f"<< Base directories are successfully created. Read below before you run the main function")
print(f"   your root ({root_dir})")
print(f"   ├── data")
print(f"   ├── input")
print(f"   ├── output")
print(f"   │   ├── tsv")
print(f"   │   ├── img")
print(f"   │   └── genus")
print(f"   └── seqlib")
