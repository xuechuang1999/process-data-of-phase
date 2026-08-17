import os
import shutil
import subprocess
import stat

# 赝势路径映射
POTCAR_PATH = {
    "Hf": os.path.expanduser("~/POTENTIAL2017/PAW_PBE/Hf_sv/POTCAR"),
    "Zr": os.path.expanduser("~/POTENTIAL2017/PAW_PBE/Zr_sv/POTCAR"),
    "O": os.path.expanduser("~/POTENTIAL2017/PAW_PBE/O/POTCAR")
}

# POSCAR 修改规则
POSCAR_RULES = {
    "a": (["Hf", "Zr", "O"], ["2", "2", "8"]),
    "b": (["Hf", "Zr", "Hf", "Zr", "O"], ["1", "1", "1", "1", "8"]),
    "c": (["Hf", "Zr", "Hf", "O"], ["1", "2", "1", "8"])
}

# run-vasp 模板路径
TEMPLATE_RUN_SCRIPT = "/gs/home/xuechuang/run-bash/run-vasp"

def modify_poscar(poscar_path, prefix):
    if not os.path.isfile(poscar_path):
        return
    lines = open(poscar_path, 'r').readlines()
    if len(lines) >= 7:
        elements, counts = POSCAR_RULES[prefix]
        lines[5] = " ".join(elements) + "\n"
        lines[6] = " ".join(counts) + "\n"
        with open(poscar_path, 'w') as f:
            f.writelines(lines)

def generate_potcar(poscar_path, potcar_path):
    # generate POTCAR according to the elements in POSCAR and the mapping in POTCAR_PATH
    if not os.path.isfile(poscar_path):
        # raise FileNotFoundError(f"POSCAR file not found: {poscar_path}")
        return
    with open(poscar_path, 'r') as f:
        lines = f.readlines()
    if len(lines) < 6:
        # raise ValueError(f"Invalid POSCAR format: {poscar_path}, the file has less than 6 lines.")
        return
    elements = lines[5].strip().split()
    with open(potcar_path, 'wb') as potcar:
        for elem in elements:
            pot_path = POTCAR_PATH.get(elem)
            if pot_path and os.path.isfile(pot_path):
                with open(pot_path, 'rb') as p:
                    shutil.copyfileobj(p, potcar)
            else:
                print(f"Warning: POTCAR for {elem} not found.")

def modify_incar(incar_path):
    """反向修改 INCAR：IALGO=48->38，IOPTCELL->#IOPTCELL，ISYM=2->-1"""
    if not os.path.isfile(incar_path):
        return
    lines = []
    with open(incar_path, 'r') as f:
        for line in f:
            if line.strip().startswith("IALGO") and "48" in line:
                lines.append("IALGO = 38\n")
            elif line.strip().startswith("IOPTCELL") and not line.strip().startswith("#"):
                lines.append("#IOPTCELL = 1 0 0 0 1 0 0 0 1\n")
            elif line.strip().startswith("ISYM") and "2" in line:
                lines.append("#ISYM   = -1\n")
            else:
                lines.append(line)
    with open(incar_path, 'w') as f:
        f.writelines(lines)

def copy_run_script(dst_dir):
    if os.path.isfile(TEMPLATE_RUN_SCRIPT):
        dst_script = os.path.join(dst_dir, "run-vasp")
        shutil.copy(TEMPLATE_RUN_SCRIPT, dst_script)
        # 添加执行权限
        st = os.stat(dst_script)
        os.chmod(dst_script, st.st_mode | stat.S_IEXEC)
    else:
        print(f"Warning: Template run-vasp script not found at {TEMPLATE_RUN_SCRIPT}")

def submit_job(directory):
    run_script = os.path.join(directory, "run-vasp")
    if os.path.isfile(run_script):
        subprocess.run(["sbatch", "run-vasp"], cwd=directory)

def process_folder(folder):
    for prefix in ["a", "b", "c"]:
        # 创建 a-/b-/c- 文件夹在原始 folder 内
        new_dir = os.path.join(folder, f"{prefix}-{folder}")
        os.makedirs(new_dir, exist_ok=True)

        for sub in ["ini", "fin"]:
            src = os.path.join(folder, sub)
            dst = os.path.join(new_dir, sub)
            if os.path.isdir(src):
                shutil.copytree(src, dst, dirs_exist_ok=True)

                # 拷贝提交脚本并加执行权限
                copy_run_script(dst)

                poscar = os.path.join(dst, "POSCAR")
                incar = os.path.join(dst, "INCAR")
                potcar = os.path.join(dst, "POTCAR")

                modify_poscar(poscar, prefix)
                generate_potcar(poscar, potcar)
                modify_incar(incar)

                # 提交任务
                submit_job(dst)

def main():
    base_dir = os.getcwd()
    for item in os.listdir(base_dir):
        path = os.path.join(base_dir, item)

        # 跳过隐藏文件夹、符号链接、a-/b-/c-开头的目录
        if not os.path.isdir(path) or item.startswith((".", "a-", "b-", "c-")) or os.path.islink(path):
            continue

        # 检查 item 内部是否已经有 a-/b-/c- 子目录
        subfolders = os.listdir(path)
        if any(name.startswith(("a-", "b-", "c-")) and os.path.isdir(os.path.join(path, name)) for name in subfolders):
            print(f"Skip {item}: already processed (contains a-/b-/c- folders)")
            continue

        # 处理文件夹
        process_folder(item)

if __name__ == "__main__":
    main()