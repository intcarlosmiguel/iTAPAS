import numpy as np
import os

Volume = 300
N_ODS = 10
# gera N_ODS volumes distintos (positivos) que somam Volume
data = np.loadtxt("./fortaleza/edges.txt")[:,:2].astype(int)
N = np.max(data)
rng = np.random.default_rng()

def gerar_volumes_distintos(total, n, max_attempts=1000):
    for _ in range(max_attempts):
        cortes = np.sort(rng.random(n - 1) * total)
        bordas = np.concatenate(([0.0], cortes, [total]))
        vols = np.diff(bordas)
        if len(np.unique(np.round(vols, 12))) == n:
            return vols
    raise RuntimeError("Não foi possível gerar volumes distintos")

def gerar_pares_od(n_pairs, n_nodes, rng):
    pairs = set()
    # tentativa aleatória inicial
    while len(pairs) < n_pairs:
        o = int(rng.integers(1, n_nodes + 1))
        d = int(rng.integers(1, n_nodes + 1))
        if o != d:
            pairs.add((o, d))
        # se for impossível obter pares únicos aleatórios (muito poucos nós),
        # preencha deterministicamente
        if len(pairs) < n_pairs and len(pairs) + (n_nodes*(n_nodes-1) - len(pairs)) <= len(pairs):
            break

    # preenchimento determinístico se necessário
    if len(pairs) < n_pairs:
        for i in range(1, n_nodes + 1):
            for j in range(1, n_nodes + 1):
                if i != j:
                    pairs.add((i, j))
                    if len(pairs) == n_pairs:
                        break
            if len(pairs) == n_pairs:
                break

    if len(pairs) < n_pairs:
        raise RuntimeError("Não foi possível gerar pares OD suficientes")
    return list(pairs)

for i in range(100):
    np.random.seed(i+42)
    od_volumes = np.array(gerar_volumes_distintos(Volume, N_ODS))
    od_pairs = np.array(gerar_pares_od(N_ODS, N, rng))
    output = "./od_outputs"
    os.makedirs(output, exist_ok=True)
    out_dir = f"./od_outputs/OD_{N_ODS}_{Volume}"
    os.makedirs(out_dir, exist_ok=True)
    OD = np.vstack((od_pairs.T, od_volumes)).T
    np.savetxt(f"{out_dir}/OD_{i}.txt", OD, fmt="%d %d %.6f")
#print(f"Volumes distintos gerados: {od_volumes}")
#print(f"{N_ODS} pares OD gerados e salvos em '{out_dir}' (arquivos OD_i.txt e OD_all.txt)")