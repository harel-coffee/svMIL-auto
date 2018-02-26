rm -rf ~/Documents/Projects/StructuralVariants/neo4j-3.3/packaging/standalone/target/neo4j-community-3.3.3-SNAPSHOT/data/databases/graph.db
sh ~/Documents/Projects/StructuralVariants/neo4j-3.3/packaging/standalone/target/neo4j-community-3.3.3-SNAPSHOT/bin/neo4j-import \
--into ~/Documents/Projects/StructuralVariants/neo4j-3.3/packaging/standalone/target/neo4j-community-3.3.3-SNAPSHOT/data/databases/graph.db \
--id-type string \
--nodes:Region Interchromosomal/regions.csv \
--relationships:interactsWith Interchromosomal/regions_regions_rel.csv \
~/Documents/Projects/StructuralVariants/neo4j-3.3/packaging/standalone/target/neo4j-community-3.3.3-SNAPSHOT/bin/neo4j restart