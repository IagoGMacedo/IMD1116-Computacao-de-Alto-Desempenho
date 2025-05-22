import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

Nx, Ny, Nz = 181, 181, 181

with open('evolucao3d.csv', 'r') as f:
    content = f.read().strip().split('\n\n')
    frames = []
    for frame in content:
        arr = np.array([float(x) for x in frame.replace('\n', ',').split(',') if x.strip() != ''])
        frames.append(arr.reshape((Nx, Ny, Nz)))

plano = 'xy'
idx = Nz // 2

def get_slice(volume, plano, idx):
    if plano == 'xy':
        return volume[:, :, idx]
    elif plano == 'xz':
        return volume[:, idx, :]
    elif plano == 'yz':
        return volume[idx, :, :]
    else:
        raise ValueError("Plano deve ser 'xy', 'xz' ou 'yz'.")

fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(get_slice(frames[0], plano, idx), origin='lower', cmap='viridis', aspect='auto')
cbar = plt.colorbar(im, ax=ax)
ax.set_title(f'Frame 0 - Corte {plano.upper()} em {idx}')

def update(frame_idx):
    data = get_slice(frames[frame_idx], plano, idx)
    im.set_data(data)
    ax.set_title(f'Frame {frame_idx} - Corte {plano.upper()} em {idx}')
    return [im]

ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=800, blit=False, repeat=True)
plt.show()