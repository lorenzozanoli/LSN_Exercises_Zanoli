import pygame
import sys

pygame.init()

width, height = 800, 400
screen = pygame.display.set_mode((width, height))
pygame.display.set_caption("Ising Model Animation")
clock = pygame.time.Clock()

with open("output/output_video/highT.dat", "r") as file:
    lines = file.readlines()

spins_list = [line.strip().split() for line in lines]

particles = []
radius = 10

for spin_sequence in spins_list:
    particles.append([])
    for i, spin in enumerate(spin_sequence):
        if spin == "1":
            color = (255, 0, 0)
        else:
            color = (0, 0, 255)
        particles[-1].append((i * 2 * radius + radius, height // 2, color))

line_index = 0

while line_index < len(particles):

    for i, particle in enumerate(particles[line_index]):
        x, y, color = particle
        pygame.draw.circle(screen, color, (x, y), radius)

    pygame.display.flip()
    clock.tick(10)

    line_index += 1

pygame.quit()