from bibpy.cfw import resolver_cfw
from bibpy.utils import carregar_rede, carregar_viagens, gerar_viagens_aleatorias
import matplotlib.pyplot as plt

rede = carregar_rede("./fortaleza/edges_fortaleza.txt")
N = rede.number_of_nodes()
viagens = gerar_viagens_aleatorias(rede, N*0.2, 1)

resultado = resolver_cfw(rede, viagens)

gaps = resultado['historico_gap']

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.semilogy(range(len(gaps)), gaps, linewidth=1.4, color='#1a73e8')
ax.axhline(1e-4, linestyle='--', color='#d93025', linewidth=0.9, label=r'$\varepsilon = 10^{-4}$')
ax.set_xlabel('Iteração $k$')
ax.set_ylabel('Gap relativo  $RE_k$')
ax.set_title('Convergência — Conjugate Frank-Wolfe (CFW)')
ax.legend()
ax.grid(True, which='both', alpha=0.3)
fig.tight_layout()
plt.savefig('gap_cfw.png', dpi=150)
plt.show()