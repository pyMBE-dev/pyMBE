from pydantic import BaseModel, Field

class PMBBaseModel(BaseModel):
    """
    Base class for all pyMBE models.

    Provides common fields and validation behavior for pyMBE templates and instances.

    Attributes:
        pmb_type (str): Fixed type identifier. Subclasses must set this to a specific type.
        name (str): Unique name of the model instance or template.

    Config:
        validate_assignment (bool): Ensures that attribute assignments are validated.
        extra (str): Forbids extra attributes not defined in the model.
    """

    pmb_type: str = Field(frozen=True)
    name: str
    
    class Config:
        validate_assignment = True
        extra = "forbid"
