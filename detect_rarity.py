import pickle
import hyperloglog

print("Cargando el catálogo R... (HLL_R.sketch)")

# 'rb' = Read Binary (Leer Binario)
with open("HLL_R.sketch", "rb") as f:
    hll_R = pickle.load(f)

print("¡Catálogo cargado!")
print(f"La cardinalidad estimada es: {len(hll_R)}")