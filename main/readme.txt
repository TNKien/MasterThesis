READ ME

Code pour calculer les réseaux de connectivité et en extraire les paramètres de connectivité.
La database fournie est composée de 100 patients aux USI (possibilité de rajouter des patients si besoin).
Pour lancer le calcul de connectivité sur les patients, run le main.

Le calcul des réseaux de connectivité se fait de la manière suivante :
1. Pre-processing des EEG 
	- re-référencement au montage Laplacien
	- filtrage dans les bandes de fréquences usuelles (delta [1-4 Hz], theta [4-8 Hz], alpha [8-12 Hz], beta [12-20 Hz], quasi-broadband [1-12 Hz] et broadband [1-20 Hz])
	- découpe des signaux en fenêtres d'une seconde, overlapping de 500 milisecondes
2. Calcul des associations entre électrodes (i.e. calcul des liens du réseau). Un ensemble de mesures sont disponibles pour ce calcul, groupées dans le code par catégories :
	- Catégorie correlation avec le maximum de cross-correlation, le lag à ce maximum, le maximum de cross-correlation corrigée et le lag à ce maximum comme mesures
	- Catégorie cohérence avec la cohérence et la partie imaginaire de la cohérence comme mesures
	- Catégorie phase locking value (PLV) avec le plv comme mesure
	- Catégorie phase avec le phase lag index (PLI) et le weighted phase lag index (wPLI) comme mesures
	- Catégorie information avec la mutual information and la weighted symbolic mutual information (wSMI) comme mesures
   Par défaut, la catégorie "phase" est spécifiée dans le main.
3. Extraction des paramètres de connectivité. Les paramètres suivants sont disponibles (vous êtes invités à compléter cette liste en fonction de vos besoins) :
	- clustering coefficient
        - characteristic path length
	- strength
	- global efficiency

Certains paramètres sont spécifiés au début du main pour le calcul de la connectivité. Vous pouvez jouer avec certains pour obtenir différentes estimations de la connectivité.
Paramètres à ne pas modifier :
	- startTime : time point du début de l'analyse
	- timeToAnalyze : nombre de time points considérés pour l'analyse
	- channels : canaux EEG
	- reference : référence d'enregistrement de l'EEG
Paramètres pouvant être modifiés :
	- fhband : bornes inférieures de chaque bande de fréquence considérée
	- flband : bornes supérieures de chaque bande de fréquence considérée
	- windowSec : taille de fenêtres
	- overlapSec : overlap entre fenêtres
	- measure : catégorie de mesure d'association (voir ci-dessus) - options : correlation(50), coherence, plv, phase, and information(200)

NB : en sortie de l'algorithme il y aura un réseau de connectivité et les paramètres associés pour chaque patient et par mesure dans la catégorie et bande de fréquence.
