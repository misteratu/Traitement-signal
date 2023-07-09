# Traitement du signal

L’objectif de ce projet était de créer un modem suivant la recommandation V21 de l’UIT (Union Internationale des Télécommunications). Ce modem aura un débit maximal de 300 bit/s. La première étape consistera à former le signal 2-FSK à transmettre à partir d’un fichier d’infor- mation binaire (modulation).
La deuxième étape est de bruiter le signal avec un bruit Gaussien pour correspondre à la réalité des transmissions.
La troisième étape est de créer un récepteur pour lire les données envoyées. Pour cela, il a été choisi de construire ce récepteur grâce à 3 méthodes différentes : La première par filtrage, la deuxième par démodulation FSK avec la synchronisation supposée idéale et la troisième par démodulation FSK avec prise en compte une erreur de synchronisation entre émetteur et récepteur.
