# help functions for mapping metadata to masterdata object types
def metadata_to_instance(metadata: dict, instance):
    """sets metadata in the object class

    Args:
        metadata (dict): Metadata of Object
        instance (_type_): Object-Class-Instance

    Returns:
        Object-Instance
    """
    props = metadata_to_masterdata(metadata, instance)
    for k, v in props.items():
        setattr(instance, k, v)
    return instance


def metadata_to_masterdata(metadata: dict, object_instance):
    """checks if avaliable metadata is in Object-Attributes

    Args:
        metadata (dict): Metadata of Object
        object_instance (_type_): Object-Class-Instance

    Returns:
         List with verified Object-Attributes
    """
    avaliable_values = object_instance._property_metadata
    data = metadata
    object_prop_list = {
        prop_name: data_value
        for prop_name, prop_obj in avaliable_values.items()
        if (data_value := data.get(prop_obj.code)) is not None
    }
    return object_prop_list
