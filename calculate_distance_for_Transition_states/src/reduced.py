#!/usr/bin/env python
import spglib
from ase.io import vasp
import numpy as np
import os


def rodrigues_rotation(r, theta):
    r = np.array(r).reshape(3, 1)
    rx, ry, rz = r[:, 0]
    M = np.array([
        [0, -rz, ry],
        [rz, 0, -rx],
        [-ry, rx, 0]
    ])
    R = np.eye(3)
    R[:3, :3] = np.cos(theta) * np.eye(3) + (1 - np.cos(theta)) * r @ r.T + np.sin(theta) * M
    return R


def align_cell(cell):
    cell = np.dot(cell, rodrigues_rotation([0, 0, 1], np.arctan(cell[0][1] / cell[0][0])))
    cell = np.dot(cell, rodrigues_rotation([0, -1, 0], np.arctan(cell[0][2] / cell[0][0])))
    cell = np.dot(cell, rodrigues_rotation([1, 0, 0], np.arctan(cell[1][2] / cell[1][1])))
    return cell


def matrix2string(rot, trans):
    M = np.sum(rot, axis=1)
    P = np.dot(np.abs(rot), np.array([1, 2, 3]))
    preffix = ''.join([chr(int(i) + 119) for i in P])
    suffix = ''.join(['1' if i == -1 else '0' for i in M])
    translation = ''.join(['0' if min(i % 1, (1 - i) % 1) < 0.1 else '1' for i in trans])
    return preffix + '-' + suffix + '-' + translation


def string2matrix(string):
    preffix, suffix, translation = string.split('-')
    M = np.diag([-1 if i == '1' else 1 for i in suffix])
    P = np.zeros((3, 3))
    for index, i in enumerate([ord(i) - 119 for i in preffix]):
        P[index, i - 1] = 1
    translation = [0.5 if i == '1' else 0 for i in translation]
    return np.dot(M, P), translation


def equal_zero_half(v, tol=1e-2):
    count = 0
    for i in range(3):
        if min(v[i] % 1, (1 - v[i]) % 1) > tol:
            if min((v[i] - 0.5) % 1, (1.5 - v[i]) % 1) <= tol:
                count += 1
            else:
                return False
    if count != 2 and count != 0:
        return False
    return True


def get_symmetry_operators(system, tol=1e-2):
    symmetry_info = spglib.get_symmetry((system.cell, system.get_scaled_positions(), system.numbers), symprec=tol)
    if not symmetry_info:
        return []
    operator_list = []
    for rot, trans in zip(symmetry_info['rotations'], symmetry_info['translations']):
        if equal_zero_half(trans, tol):
            operator_list.append(matrix2string(rot, trans))
        else:
            raise ValueError('Translation vector is not zero or half of the cell vector.')
    operator_list.sort()
    return operator_list


def reduced(system):
    symbol_list = []
    for i in ['xyz', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx']:
        for j in ['000', '001', '010', '011', '100', '101', '110', '111']:
            for k in ['000', '011', '101', '110']:
                symbol_list.append(i + '-' + j + '-' + k)
    states = []
    rot_symbol = get_symmetry_operators(system)
    states.append(rot_symbol)
    for operator in symbol_list:
        for group in states:
            if operator in group:
                break
        else:
            new_state = []
            for i in rot_symbol:
                rot1, trans1 = string2matrix(operator)
                rot2, trans2 = string2matrix(i)
                trans2 = np.dot(rot1, np.array(trans2))
                new_state.append(matrix2string(np.dot(rot1, rot2), [i + j for i, j in zip(trans1, trans2)]))
            states.append(new_state)
    return states


def process(system):
    state_list = reduced(system)
    with open('states', 'w') as f:
        for i in state_list:
            for j in i:
                f.write(j + ' ')
            f.write('\n')
    with open('states1', 'w') as f:
        for i in state_list:
            for j in i:
                ls = []
                for k in range(3):
                    string = ''
                    if j[k + 4] == '1':
                        string+='-'
                    string+=j[k]
                    if j[k + 8] =='1':
                        string+='+1/2'
                    ls.append(string)
                f.write('(' + ','.join(ls) + ') ')
            f.write('\n')
    os.makedirs('reduced', exist_ok=True)
    for state in state_list:
        state_dir = os.path.join('reduced', state[0])
        os.makedirs(state_dir, exist_ok=True)
        new_cell = system.cell[:, [ord(i) - 120 for i in state[0][:3]]][[ord(i) - 120 for i in state[0][:3]], :]
        for i in range(3):
            if state[0][i + 4] == '1':
                new_cell[i, :] = -new_cell[i, :]
                new_cell[:, i] = -new_cell[:, i]
        new_cell = align_cell(new_cell)
        new_positions = np.transpose([1 - system.get_scaled_positions()[:, ord(state[0][i]) - 120] if state[0][i + 4] == '1' else system.get_scaled_positions()[:, ord(state[0][i]) - 120] for i in range(3)])
        new_positions = new_positions + np.array([float(i) * 0.5 for i in state[0][8:11]])
        new_positions = new_positions % 1
        new_positions = sort_positions(new_positions)
        new_system = system.copy()
        new_system.set_cell(new_cell)
        new_system.set_scaled_positions(new_positions)
        new_system.write(os.path.join(state_dir, 'POSCAR'), vasp5=True, direct=True)
        os.system('cp INCAR KPOINTS POTCAR run-vasp '+state_dir)
        os.system('cd '+state_dir+'; sbatch run-vasp')


def sort_positions(positions):
    hf = positions[:4]
    hf1 = sorted(hf, key=lambda x: dis(x, [0.0, 0.0, 0.0]))[0]
    hf2 = sorted(hf, key=lambda x: dis(x, [0.0, 0.5, 0.5]))[0]
    hf3 = sorted(hf, key=lambda x: dis(x, [0.5, 0.0, 0.5]))[0]
    hf4 = sorted(hf, key=lambda x: dis(x, [0.5, 0.5, 0.0]))[0]
    o = sorted(positions[4:], key=lambda x: [i <= 0.5 for i in x])
    return [hf1, hf2, hf3, hf4] + o


def dis(coord1, coord2):
    a = 0
    for i in range(3):
        a += min((coord1[i] - coord2[i])**2, ((coord1[i] - coord2[i]) % 1 - 1)**2, ((coord1[i] - coord2[i]) % 1)**2)
    return a**0.5


if __name__ == '__main__':
    process(vasp.read_vasp('POSCAR'))