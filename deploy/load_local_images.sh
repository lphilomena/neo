
IMAGE_DIR="${NEOAG_CONTAINER_IMAGE_DIR:-/path/to/neodata4git/work/container_images}"

docker load -i "${IMAGE_DIR}/neoag-hla-la_ubuntu22.04.tar"
docker load -i "${IMAGE_DIR}/neoag-easyfuse_ubuntu22.04.tar"
#docker load -i "${IMAGE_DIR}/neoag-netmhcpan_4.2c-ubuntu22.04.tar"
docker load -i "${IMAGE_DIR}/neoag-netmhcstabpan_1.0-ubuntu22.04.tar"
docker load -i "${IMAGE_DIR}/neoag-purple-suite_ubuntu22.04.tar"
docker load -i "${IMAGE_DIR}/neoag-spechla_ubuntu22.04.tar"
