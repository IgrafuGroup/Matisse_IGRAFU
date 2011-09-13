#IGRAFU

A IGRAFU (Interface Gráfica para Reconstruções de Arvores Filogenéticas da UESC) visa agrupar programas open-source gratuitos que realizam a reconstrução
de arvores filogenéticas e/ou dão suporte a esta. Além disso a IGRAFU inova com uma GUI de usuário de fácil aprendizado, porém robusta e diversificada.

Objetivo da IGRAFU é prover, ao cientista que deseja realizar a reconstrução de arvores filogenéticas, um ambiente unificado com uma curva de aprendizado menor, mas que ofereça todos os recursos dos atuais programas em linha de comando. Por fim o IGRAFU é software livre e gratuito diferente de ferramentas similares.

##Matisse_IGRAFU

Esta versão aqui apresentada - Matisse_Igrafu- é uma implementação realizada através do NetBeans GUI Builder (apelidado Matisse). 
Nesta presente distribuição encontra-se presente funcionalmente somente o método de inferência Bayesiana através de uma GUI para o MrBayes-3.1.2.
Se utiliza também do FigTree para visualização das arvores geradas.

#Instalação

A instalação requer somente a compilação do mrbayes-3.1.2. 
Existe contudo um script de instalação que automatiza esta tarefa. Segue os passo para executá-lo:

	% chmod +x Install
	% ./Install

#Execução

Para executar o Matisse_IGRAFU basta executar o script IGRAFU:

	% ./IGRAFU &


