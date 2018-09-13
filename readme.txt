Även om jag har varit i princip alla filer och pillat på gammal Fortran 66 och en del
annat så ligger all det nya med egenvärdeslösaren i mappen rscf_test och filen maneig.f.
JUst nu är allt väldigt enkelt kodat med en parameter solver som sätts i början av filen.
Om solver sätts till 1 så kommer FEAST att användas. För varje block kommer först ett anrop
till subrutinen iniest3 att göras. Iniest3 är en variant av den vanliga rutinen för att
uppskatta egenvärdena, med det nya att den bara svarar med en grov uppskattning av det lägst 
och högst liggande egenvärdet. Efter anropet till iniest3 frigörs minnet från alla tillfälliga
arbetsvariabler.

FEAST arbetar till skillnad från DVDSON inom ett givet intervall. Det visade sig att iniest
oftast ganska grovt överskattar intervallet vilket ger ett onödigt stort antal egenvärden
som ska hittas. Därför inleds lösningen med att FEAST först anropas med uppdraget att ge
en uppskattning av antalet egenvärden, M, i det givna intervallet. Sedan minskas intervallet
uppifrån iterativt tills M närmar sig det efterfrågade antalet egenvärden. Detta görs
med begränsad numerisk upplösning (hela matrisen men färre noder för integration) så det
går relativt fort. Den interna dimensioneringsparametern M0 sätts till ca 1.5*M och
variablerna dimensioneras därefter med Fortran90-alloc. FEAST har en intern konvergensparameter
som just nu är satt till 12 (spåret mindre än 10^-12) vilket ger fler konvergerade decimaler
än DVDSON så där går det nog att spara lite tid. Avslutningsvis samlas de konvergerade
egenvärdena och -vektorerna ihop som tidigare. För att kunna göra den preliminära uppskattningen
krävs FEAST 3.0 eller senare. Intels MKL använder FEAST 2.1 så FEAST måste kompileras
från grunden. Den bifogade koden är enbart ändrad på ett ställe märkt OOC där en parameter
sätts explicit, se nedan.

FEAST har flera möjliga moder för att lösa egenvärdesproblemet och jag har valt den lättaste
utan omvänd kommunikation som använder PARDISO för diagonalisering, men det går även att
använda omvänd kommunikation och en valfri diagonaliserare. Jag har valt att använda MKLs
version av PARDISO för den är mycket snabb och bra parallelliserad. Notera att MKL fram
till och med 11.1.4 har en bug som gör att PARDISO kraschar för stora matriser. Jag använder
11.3 med god framgång. Det går utmärkt att kompilera PARDISO och de stödbibliotek som behövs
till PARDISO från källkod och alla delar går att få tag i gratis under GPL och motsvarande
licenser på nätet.

MKLs PARDISO kan med parametern som jag kallat OOC i en kommentar i FEAST-koden fås att
antingen alltid köra i minnet eller alltid mot disk eller som en hybrid som flyttar ut på
disk om den interna minnesåtgången överstiger värdet satt av shell-variabeln
MKL_PARDISO_OOC_MAX_CORE_SIZE som sätter gränsen i MB. Defaultvärdet är 2000.

Eftersom jag mest har varit intresserad av att testa konceptet har jag bara utnyttjat
OpenMP-parallelliseringen av PARDISO men det går att parallellisera på flera nivåer. Det
går att kompilera en egen version av PARDISO (med stödbibliotek) som använder både
MPI och OpenMP för att köra över flera separata datorer, det går att kompilera en MPI-
parallell version av FEAST och det går att i GRASP dela upp energiintervallet i flera
delintervall som sedan löses som separata delproblem. Idag utnyttjar jag ingen
parallellism i GRASP men det går såklart också att göra med lite kodande. Jag har en del
idéer som jag inte har hunnit testa.

För att slippa trassla med att länka in MKL med allt som behövs med hjälp av gfortran så
har jag komplierat allt med ifort och resultatet blir detsamma som för gfortran på
"den vanliga" delen av koden. Med ifort 16.0.0 så går det förvisso ca 30% snabbare än med
gfortran, men det kan jag stå ut med :-)
